"""
analysis_common.py -- shared helpers for the SQ2 spatial-GRN analysis.

Loads the cached spatial-link matrices (built by 13_spatial_link.py) and provides
the building blocks every test needs:
  * per-subtype row selection
  * within-section pair blocks (spatial-distance + GRN-dissimilarity matrices)
  * within-section kNN spatial-weights matrix (block-diagonal by section)
  * Moran's I

KEY METHOD CONSTRAINTS (forced by the data / established decisions):
  * Coordinates are argmax (x, y). Each MERFISH *section* is its own coordinate
    frame (not co-registered), so spatial distance is ONLY ever computed
    WITHIN a section. There is no depth-from-pia column available, so the
    spatial axis is within-section Euclidean distance.
  * Everything is run PER SUBTYPE (cell-type is a confound otherwise).
  * Expression-space versions of each test are computed as the circularity
    baseline (Tangram places cells by expression similarity, so a GRN-vs-space
    signal must be shown to exceed the raw expression-vs-space signal).
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.neighbors import NearestNeighbors

DFLAM    = Path("/tudelft.net/staff-umbrella/ScReNI/dflam")
CACHE    = DFLAM / "data/spatial_link"
ANALYSIS = DFLAM / "analysis"

SUBTYPES = ["Astrocyte", "L2/3 IT", "L4 IT", "Oligodendrocyte"]
SLUG = {"Astrocyte": "astrocyte", "L2/3 IT": "l23_it",
        "L4 IT": "l4_it", "Oligodendrocyte": "oligo"}

# analysis hyper-parameters (single source of truth)
MIN_CELLS_PER_SECTION = 30      # sections smaller than this are skipped
MAX_CELLS_PER_SECTION = 1500    # subsample cap to keep pairwise O(n^2) cheap
KNN_K                 = 6       # within-section spatial neighbours for W / Moran
N_PERM                = 999     # permutations for p-values
N_PCS                 = 20      # GRN/expression PCA components
SEED                  = 42


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------
def load_cache(with_expr=False):
    """Return (X_grn[csr cells x edges], cell_ids, meta[aligned], expr_or_None)."""
    X = sp.load_npz(CACHE / "spatial_grn_matrix.npz").tocsr()
    idx = np.load(CACHE / "spatial_grn_index.npz", allow_pickle=True)
    cell_ids = idx["cell_ids"].astype(str)
    meta = pd.read_parquet(CACHE / "spatial_meta.parquet")
    meta = meta.loc[cell_ids].copy()          # align to matrix row order
    meta["section"] = meta["section"].astype(str)
    expr = None
    if with_expr:
        ez = np.load(CACHE / "spatial_expr.npz", allow_pickle=True)
        e_ids = ez["cell_ids"].astype(str)
        E = np.asarray(ez["expr"], dtype=np.float32)
        # align expression rows to the GRN cell order
        pos = pd.Series(np.arange(len(e_ids)), index=e_ids)
        E = E[pos.loc[cell_ids].values]
        expr = E
    if X.shape[0] != len(cell_ids):
        raise ValueError(f"matrix rows {X.shape[0]} != cell_ids {len(cell_ids)}")
    return X, cell_ids, meta, expr


def subtype_rows(meta, subtype):
    """Row indices (into the full matrix) for one subtype."""
    return np.flatnonzero(meta["subtype"].values == subtype)


# --------------------------------------------------------------------------
# within-section pair blocks (variogram + Mantel)
# --------------------------------------------------------------------------
def _cosine_dissim(M):
    """1 - cosine similarity, as a dense (n x n) float32 matrix."""
    S = cosine_similarity(M)
    np.clip(S, -1.0, 1.0, out=S)
    D = (1.0 - S).astype(np.float32)
    np.fill_diagonal(D, 0.0)
    return D


def section_blocks(feature_matrix, meta, rows, rng,
                   min_cells=MIN_CELLS_PER_SECTION,
                   max_cells=MAX_CELLS_PER_SECTION):
    """
    Split one subtype's cells into per-section blocks and, for each section,
    build the spatial-distance and feature-dissimilarity square matrices.

    feature_matrix : either the sparse GRN matrix (rows = all cells) OR a dense
                     (n_cells x d) array (e.g. expression). Indexed by `rows`.
    Returns list of dicts: {section, n, rows(global), S(spatial), D(dissim)}.
    """
    coords = meta[["x", "y"]].values
    sect = meta["section"].values
    blocks = []
    for s in pd.unique(sect[rows]):
        sel = rows[sect[rows] == s]
        if len(sel) < min_cells:
            continue
        if len(sel) > max_cells:
            sel = rng.choice(sel, size=max_cells, replace=False)
            sel.sort()
        xy = coords[sel]
        # spatial distance (within this section's own frame)
        diff = xy[:, None, :] - xy[None, :, :]
        S = np.sqrt((diff ** 2).sum(-1)).astype(np.float32)
        # feature dissimilarity
        if sp.issparse(feature_matrix):
            D = _cosine_dissim(feature_matrix[sel])
        else:
            D = _cosine_dissim(feature_matrix[sel])
        blocks.append({"section": s, "n": len(sel), "rows": sel, "S": S, "D": D})
    return blocks


def pooled_pairs(blocks):
    """Concatenate the upper-triangle pairs across all section blocks."""
    sd, dd = [], []
    for b in blocks:
        iu = np.triu_indices(b["n"], k=1)
        sd.append(b["S"][iu])
        dd.append(b["D"][iu])
    return np.concatenate(sd), np.concatenate(dd)


# --------------------------------------------------------------------------
# within-section spatial weights + Moran's I
# --------------------------------------------------------------------------
def within_section_knn_W(meta, rows, k=KNN_K):
    """
    Row-standardised kNN spatial-weights matrix (sparse, n x n where
    n = len(rows)), neighbours drawn ONLY from the same section.
    """
    coords = meta[["x", "y"]].values[rows]
    sect = meta["section"].values[rows]
    n = len(rows)
    rptr, cidx = [], []
    for s in pd.unique(sect):
        loc = np.flatnonzero(sect == s)
        if len(loc) <= 1:
            continue
        kk = min(k, len(loc) - 1)
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(coords[loc])
        _, nbr = nn.kneighbors(coords[loc])
        for a in range(len(loc)):
            for b in nbr[a, 1:]:            # skip self (first neighbour)
                rptr.append(loc[a]); cidx.append(loc[b])
    data = np.ones(len(rptr), dtype=np.float64)
    W = sp.csr_matrix((data, (rptr, cidx)), shape=(n, n))
    rs = np.asarray(W.sum(1)).ravel()
    rs[rs == 0] = 1.0
    W = sp.diags(1.0 / rs) @ W            # row-standardise
    return W


def morans_I(z, W):
    """Global Moran's I for a single (already mean-centred) feature vector z."""
    n = len(z)
    Wsum = W.sum()
    num = float(z @ (W @ z))
    den = float(z @ z)
    if den == 0 or Wsum == 0:
        return np.nan
    return (n / Wsum) * (num / den)


def morans_I_perm(feature, W, rng, n_perm=N_PERM):
    """Moran's I + two-sided permutation p-value for one feature."""
    z = feature - feature.mean()
    I_obs = morans_I(z, W)
    n = len(z)
    perm = np.empty(n_perm)
    for p in range(n_perm):
        zp = z[rng.permutation(n)]
        perm[p] = morans_I(zp, W)
    # two-sided: how often |perm - E| >= |obs - E|, E = mean of perm null
    E = perm.mean()
    extreme = np.sum(np.abs(perm - E) >= np.abs(I_obs - E))
    p_val = (extreme + 1) / (n_perm + 1)
    return I_obs, E, p_val


# --------------------------------------------------------------------------
# stratified (within-section) Mantel
# --------------------------------------------------------------------------
def _pearson(x, y):
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x @ x) * (y @ y))
    return float((x @ y) / d) if d > 0 else np.nan


def stratified_mantel(blocks, rng, n_perm=N_PERM):
    """
    Mantel correlation between spatial distance and feature dissimilarity,
    pooling pairs across sections. Significance via permutation of cell labels
    *within each section* (restricted permutation respecting the block frames).

    Returns dict with pearson r, spearman r, permutation p, n_pairs, n_sections.
    """
    from scipy.stats import rankdata
    sd, dd = pooled_pairs(blocks)
    r_obs = _pearson(sd, dd)
    r_spear = _pearson(rankdata(sd), rankdata(dd))

    # observed pooled cross-statistic is monotone in r given fixed sd; we
    # permute within section and recompute the pooled pearson each time.
    perm_r = np.empty(n_perm)
    for p in range(n_perm):
        dd_perm = []
        for b in blocks:
            order = rng.permutation(b["n"])
            Dp = b["D"][np.ix_(order, order)]
            iu = np.triu_indices(b["n"], k=1)
            dd_perm.append(Dp[iu])
        dd_perm = np.concatenate(dd_perm)
        perm_r[p] = _pearson(sd, dd_perm)
    extreme = np.sum(np.abs(perm_r) >= np.abs(r_obs))
    p_val = (extreme + 1) / (n_perm + 1)
    return {"r_pearson": r_obs, "r_spearman": r_spear, "p_perm": p_val,
            "n_pairs": int(len(sd)), "n_sections": len(blocks)}


def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)
    return Path(p)
