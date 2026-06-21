"""
SQ2 / variogram robustness -- weighted Jaccard (Ruzicka) instead of cosine.

Confirms the (flat) gradient result is not an artifact of the cosine metric.
Weighted Jaccard similarity = sum(min(a,b)) / sum(max(a,b)) over edges; we plot
dissimilarity = 1 - that, vs within-section spatial distance.

For non-negative vectors, sum(min)/sum(max) is computed efficiently via sparse
element-wise minimum on sampled within-section pairs (all-pairs is O(n^2); a large
random sample of pairs is ample for a variogram).

Outputs -> analysis/SQ2/variogram/
  variogram_wjaccard_<slug>.csv   per-bin dist_center, mean dissim, sem, n_pairs
  variogram_wjaccard_all.png      all subtypes
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "variogram")
N_BINS = 20
MAX_PAIRS_PER_SECTION = 60000     # sampled within-section pairs per section


def wjaccard_dissim(M, ii, jj, chunk=2000):
    """1 - weighted Jaccard for row pairs (ii, jj) of sparse non-negative M.

    Processed in chunks: A.minimum(B) on a full block of tens of thousands of
    rows x 246,890 edges blows up memory (billions of nonzeros), so we bound it.
    """
    out = np.empty(len(ii), dtype=np.float64)
    for start in range(0, len(ii), chunk):
        sl = slice(start, start + chunk)
        A = M[ii[sl]]
        B = M[jj[sl]]
        smin = np.asarray(A.minimum(B).sum(axis=1)).ravel()
        sa = np.asarray(A.sum(axis=1)).ravel()
        sb = np.asarray(B.sum(axis=1)).ravel()
        smax = sa + sb - smin
        out[sl] = 1.0 - np.where(smax > 0, smin / smax, 1.0)
    return out


def binned(sd, dd, edges):
    which = np.digitize(sd, edges) - 1
    c, m, e, n = [], [], [], []
    for b in range(len(edges) - 1):
        mm = which == b
        if mm.sum() < 10:
            continue
        d = dd[mm]
        c.append(0.5 * (edges[b] + edges[b + 1]))
        m.append(d.mean()); e.append(d.std(ddof=1) / np.sqrt(len(d))); n.append(int(mm.sum()))
    return np.array(c), np.array(m), np.array(e), np.array(n)


def main():
    rng = np.random.default_rng(ac.SEED)
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    # weighted Jaccard assumes non-negative weights
    Xnn = X.copy()
    Xnn.data = np.clip(Xnn.data, 0, None)
    coords = meta[["x", "y"]].values
    sect = meta["section"].values
    summary = {}

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        sd_all, dd_all = [], []
        nsec = 0
        for s in pd.unique(sect[rows]):
            sel = rows[sect[rows] == s]
            if len(sel) < ac.MIN_CELLS_PER_SECTION:
                continue
            nsec += 1
            n = len(sel)
            n_all = n * (n - 1) // 2
            npairs = int(min(n_all, MAX_PAIRS_PER_SECTION))
            ii = rng.integers(0, n, size=npairs * 2)
            jj = rng.integers(0, n, size=npairs * 2)
            keep = ii != jj
            ii, jj = ii[keep][:npairs], jj[keep][:npairs]
            d = wjaccard_dissim(Xnn, sel[ii], sel[jj])
            xy = coords[sel]
            dist = np.sqrt(((xy[ii] - xy[jj]) ** 2).sum(axis=1))
            sd_all.append(dist); dd_all.append(d)
        sd = np.concatenate(sd_all); dd = np.concatenate(dd_all)
        edges = np.unique(np.quantile(sd, np.linspace(0, 1, N_BINS + 1)))
        c, m, e, ncnt = binned(sd, dd, edges)
        pd.DataFrame({"dist_center": c, "wjac_mean_dissim": m,
                      "wjac_sem": e, "n_pairs": ncnt}).to_csv(
            OUT / f"variogram_wjaccard_{ac.SLUG[st]}.csv", index=False)
        summary[st] = (c, m)
        print(f"[{st}] {nsec} sections | {len(sd):,} sampled pairs | "
              f"wJaccard dissim {m.min():.4f} -> {m.max():.4f} (range {m.max()-m.min():.4f})")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    for i, (st, (c, m)) in enumerate(summary.items()):
        ax.plot(c, m, "-o", ms=3, lw=1.5, color=f"C{i}", label=st)
    ax.set_xlabel("within-section spatial distance (µm)")
    ax.set_ylabel("mean GRN dissimilarity (1 - weighted Jaccard)")
    ax.set_title("SQ2 variogram (weighted Jaccard) — all subtypes")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "variogram_wjaccard_all.png", dpi=150)
    plt.close(fig)
    print("variogram_wjaccard done ->", OUT)


if __name__ == "__main__":
    main()
