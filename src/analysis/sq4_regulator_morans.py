"""
SQ4 / regulator_morans -- Stage 1 wide screen (the one new statistic).

For every regulator R (= a gene that appears as edge_src, i.e. a TF-candidate whose
row of the per-cell 500x500 network is its regulon), ask whether R's regulon is
SPATIALLY ORGANISED across cells, capturing magnitude AND shape jointly in one number.

Method (frozen in SQ4_plan.md S2):
  * per cell i, R's regulon vector r_i = X[i, cols(edge_src==R)]  (signed weights).
  * GMM / signed generalized-Jaccard similarity between two regulon vectors:
        inter = sum_t min(|a_t|,|b_t|) * sign(a_t) * sign(b_t)
        union = sum_t max(|a_t|,|b_t|)
        K = inter/union                  (union==0  =>  K=1  : off ~ off)
    -> bounded [-1,1], K(a,a)=1, K(a,-a)=-1, magnitude- and shape-sensitive.
  * kernel (multivariate) Moran's I, PER SECTION (each section its own frame, >=30 cells):
        K~ = H K H              (double-centre within section, H = I - 11'/n)
        num_s = sum_ij W_ij K~_ij        (W = within-section 6-NN, row-standardised)
        den_s = trace(K~)
    pooled  I_R = (sum_s num_s) / (sum_s den_s)        (variance-weighted over sections).
  * significance: within-section label permutation -> reindex the precomputed K~
    (positions/W fixed; den_s = trace is permutation-invariant). 999 perms.
        two-sided p = (1 + #{|I_perm| >= |I_obs|}) / (1 + 999).
  * NO activity filter. The only exclusion is the numeric guard: a regulator with no
    variation in a subtype (sum_s den_s == 0) is reported as NaN (not 0) and dropped
    from FDR. Every row carries n_active.

Runs as a SLURM ARRAY over regulators (strided shards). Each task writes one shard CSV;
a final --aggregate pass concatenates all shards and adds BH-FDR (q) per subtype.

Outputs -> analysis/SQ4/
  shards/sq4_regulator_morans_shard<NNN>.csv   (per array task)
  sq4_regulator_morans.csv                     (aggregated + q, after --aggregate)
"""
import argparse
import numpy as np
import pandas as pd

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ4")
SHARD_DIR = ac.ensure_dir(OUT / "shards")


# --------------------------------------------------------------------------
# GMM / signed generalized-Jaccard kernel
# --------------------------------------------------------------------------
def gmm_kernel(M, block=48):
    """
    Pairwise signed generalized-Jaccard (GMM) similarity for the rows of M.

    M : (n, d) float, signed regulon vectors (one row per cell).
    Returns dense symmetric (n, n) float64 K, K[i,j] in [-1,1], diag = 1.
    Computed in row-blocks to cap the (block, n, d) broadcast memory.
    """
    n, d = M.shape
    absM = np.abs(M).astype(np.float64)
    sgnM = np.sign(M).astype(np.float64)
    K = np.ones((n, n), dtype=np.float64)
    for i0 in range(0, n, block):
        i1 = min(i0 + block, n)
        aba = absM[i0:i1, None, :]                      # (b,1,d)
        sga = sgnM[i0:i1, None, :]                       # (b,1,d)
        mx = np.maximum(aba, absM[None, :, :])           # (b,n,d)
        union = mx.sum(-1)                               # (b,n)
        del mx
        mn = np.minimum(aba, absM[None, :, :])           # (b,n,d)
        inter = (mn * (sga * sgnM[None, :, :])).sum(-1)  # (b,n)
        del mn
        with np.errstate(invalid="ignore", divide="ignore"):
            blk = np.where(union > 0.0, inter / union, 1.0)
        K[i0:i1] = blk
    return K


def _unit_tests():
    """Assert the kernel's defining properties before any heavy work (S7 guard 1)."""
    a = np.array([[2.0, -1.0, 0.0, 3.0]])
    # K(a,a)=1
    assert np.isclose(gmm_kernel(np.vstack([a, a]))[0, 1], 1.0)
    # K(a,-a)=-1
    assert np.isclose(gmm_kernel(np.vstack([a, -a]))[0, 1], -1.0)
    # K([2,2,2],[1,1,1]) = 3/6 = 0.5  (magnitude-sensitive)
    m = np.array([[2.0, 2.0, 2.0], [1.0, 1.0, 1.0]])
    assert np.isclose(gmm_kernel(m)[0, 1], 0.5)
    # K(0,0)=1  (off ~ off)
    z = np.zeros((2, 3))
    assert np.isclose(gmm_kernel(z)[0, 1], 1.0)
    # K(0,a)=0  (silent vs active)
    za = np.array([[0.0, 0.0, 0.0], [1.0, -2.0, 3.0]])
    assert np.isclose(gmm_kernel(za)[0, 1], 0.0)
    # symmetry + bounded + unit diagonal on a random signed matrix
    rng = np.random.default_rng(0)
    R = rng.standard_normal((20, 7)) * (rng.random((20, 7)) > 0.5)  # sparse signed
    Kr = gmm_kernel(R)
    assert np.allclose(Kr, Kr.T)
    assert np.all(Kr <= 1.0 + 1e-9) and np.all(Kr >= -1.0 - 1e-9)
    assert np.allclose(np.diag(Kr), 1.0)
    print("[unit tests] GMM kernel properties OK", flush=True)


# --------------------------------------------------------------------------
# per (regulator, subtype) kernel-Moran
# --------------------------------------------------------------------------
def regulator_subtype(X, meta, rows, cols, rng, n_perm=ac.N_PERM):
    """
    Pooled GMM kernel-Moran's I + permutation p for one regulator in one subtype.
    rows : global row indices of this subtype's cells.
    cols : matrix columns of this regulator (edge_src == R).
    Returns a dict of summary stats (I_R NaN when no variation).
    """
    sect = meta["section"].values
    sub_sect = sect[rows]

    sum_num = 0.0
    sum_den = 0.0
    perm_num_total = np.zeros(n_perm, dtype=np.float64)
    # equal-section (section-balanced) accumulators -- sensitivity check #1:
    # each section gets ONE equal vote, so a single large/variable section cannot
    # carry the result. Compared against the variance-weighted I_R in aggregate().
    eqsec_Is = []
    perm_Is_total = np.zeros(n_perm, dtype=np.float64)
    n_eqsec = 0
    DEN_FLOOR = 1e-12                  # skip sections too flat to give a stable I_s
    n_active = 0
    n_cells_used = 0
    used_sections = 0

    for s in pd.unique(sub_sect):
        sel = rows[sub_sect == s]
        if len(sel) < ac.MIN_CELLS_PER_SECTION:
            continue
        n_s = len(sel)

        # regulon submatrix for this section (cells x active targets), dense + signed
        M = X[sel][:, cols].toarray().astype(np.float64)
        n_active += int((np.abs(M).sum(1) > 0).sum())
        n_cells_used += n_s

        K = gmm_kernel(M)
        # double-centre within the section: K~ = H K H  (K symmetric => row mean = col mean)
        rm = K.mean(0)
        gm = K.mean()
        Ktil = K - rm[None, :] - rm[:, None] + gm
        den_s = float(np.trace(Ktil))

        # within-section 6-NN row-standardised spatial weights (local order == sel order)
        Ws = ac.within_section_knn_W(meta, sel).tocoo()
        wr, wc, wd = Ws.row, Ws.col, Ws.data
        num_s = float((Ktil[wr, wc] * wd).sum())

        sum_num += num_s
        sum_den += den_s
        used_sections += 1

        # within-section permutation: reindex K~ (positions/W fixed, den invariant)
        perms = np.argsort(rng.random((n_perm, n_s)), axis=1)  # (P, n_s) permutations
        mapped = Ktil[perms[:, wr], perms[:, wc]]              # (P, nnz)
        pn_s = mapped @ wd                                     # (P,) section perm numerator
        perm_num_total += pn_s

        # equal-section contribution (one vote per measurable section)
        if den_s > DEN_FLOOR:
            eqsec_Is.append(num_s / den_s)
            perm_Is_total += pn_s / den_s
            n_eqsec += 1

    if sum_den == 0.0:
        return {"I_R": np.nan, "E_perm": np.nan, "p_perm": np.nan,
                "I_R_eqsec": np.nan, "p_eqsec": np.nan, "n_eqsec": n_eqsec,
                "n_active": n_active, "n_cells_used": n_cells_used,
                "n_sections": used_sections,
                "I_s_min": np.nan, "I_s_median": np.nan, "I_s_max": np.nan}

    I_obs = sum_num / sum_den
    I_perm = perm_num_total / sum_den
    E = float(I_perm.mean())
    extreme = int(np.sum(np.abs(I_perm) >= np.abs(I_obs)))
    p = (extreme + 1) / (n_perm + 1)

    # equal-section pooled statistic + its own within-section permutation p
    if n_eqsec > 0:
        I_eqsec = float(np.mean(eqsec_Is))
        I_perm_eq = perm_Is_total / n_eqsec
        ex_eq = int(np.sum(np.abs(I_perm_eq) >= np.abs(I_eqsec)))
        p_eqsec = (ex_eq + 1) / (n_perm + 1)
    else:
        I_eqsec, p_eqsec = np.nan, np.nan

    Is = np.array(eqsec_Is) if eqsec_Is else np.array([np.nan])
    return {"I_R": float(I_obs), "E_perm": E, "p_perm": float(p),
            "I_R_eqsec": I_eqsec, "p_eqsec": p_eqsec, "n_eqsec": n_eqsec,
            "n_active": n_active, "n_cells_used": n_cells_used,
            "n_sections": used_sections,
            "I_s_min": float(np.nanmin(Is)), "I_s_median": float(np.nanmedian(Is)),
            "I_s_max": float(np.nanmax(Is))}


# --------------------------------------------------------------------------
# BH-FDR (per subtype, over testable regulators)
# --------------------------------------------------------------------------
def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return np.array([])
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0.0, 1.0)
    return out


def aggregate():
    """Concatenate shard CSVs, add BH-FDR q per subtype for BOTH poolings, write the
    master CSV, and report variance-weighted vs equal-section concordance (sensitivity #1)."""
    from scipy.stats import spearmanr
    shards = sorted(SHARD_DIR.glob("sq4_regulator_morans_shard*.csv"))
    if not shards:
        raise SystemExit("no shard CSVs found to aggregate")
    df = pd.concat([pd.read_csv(f) for f in shards], ignore_index=True)
    df = df.drop_duplicates(subset=["subtype", "regulator"]).reset_index(drop=True)
    df["q_fdr"] = np.nan
    df["q_eqsec"] = np.nan
    for st in df["subtype"].unique():
        m = (df.subtype == st) & df.p_perm.notna()
        df.loc[m, "q_fdr"] = bh_fdr(df.loc[m, "p_perm"].values)
        me = (df.subtype == st) & df.p_eqsec.notna()
        df.loc[me, "q_eqsec"] = bh_fdr(df.loc[me, "p_eqsec"].values)
    df = df.sort_values(["subtype", "I_R"], ascending=[True, False])
    out = OUT / "sq4_regulator_morans.csv"
    df.to_csv(out, index=False)
    n_reg = df["regulator"].nunique()
    print(f"[aggregate] {len(shards)} shards -> {len(df)} rows, "
          f"{n_reg} regulators -> {out}", flush=True)

    # variance-weighted vs equal-section concordance (sensitivity #1):
    # if the two rankings agree, the spatial signal is not driven by one big section.
    TOPK = 20
    concord = []
    for st in df["subtype"].unique():
        g = df[(df.subtype == st) & df.I_R.notna() & df.I_R_eqsec.notna()]
        testable = int(g.shape[0])
        rho = float(spearmanr(g.I_R, g.I_R_eqsec).correlation) if testable > 2 else np.nan
        top_w = set(g.sort_values("I_R", ascending=False).head(TOPK).regulator)
        top_e = set(g.sort_values("I_R_eqsec", ascending=False).head(TOPK).regulator)
        overlap = len(top_w & top_e)
        concord.append({"subtype": st, "testable": testable,
                        "spearman_weighted_vs_equal": rho,
                        "topk": min(TOPK, testable), "topk_overlap": overlap,
                        "q_fdr_lt05": int((g.q_fdr < 0.05).sum()),
                        "q_eqsec_lt05": int((g.q_eqsec < 0.05).sum())})
        print(f"  [{st}] testable={testable} | weighted q<.05={int((g.q_fdr < 0.05).sum())} "
              f"eq-sec q<.05={int((g.q_eqsec < 0.05).sum())} | "
              f"Spearman(weighted,equal)={rho:.3f} | "
              f"top{TOPK} overlap={overlap}/{min(TOPK, testable)}")
    pd.DataFrame(concord).to_csv(OUT / "sq4_pooling_concordance.csv", index=False)


# --------------------------------------------------------------------------
# main (one array shard)
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    ap.add_argument("--aggregate", action="store_true",
                    help="concatenate shard CSVs + add BH-FDR, then exit")
    args = ap.parse_args()

    _unit_tests()

    if args.aggregate:
        aggregate()
        return

    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
    edge_src = idx["edge_src"].astype(str)
    if edge_src.shape[0] != X.shape[1]:
        raise ValueError(f"edge_src {edge_src.shape[0]} != n_cols {X.shape[1]}")

    regulators = np.unique(edge_src)            # sorted, deterministic (expect 500)
    my_regs = regulators[args.shard::args.n_shards]
    print(f"[shard {args.shard}/{args.n_shards}] {len(my_regs)} regulators "
          f"of {len(regulators)} | matrix {X.shape}", flush=True)

    # pre-slice subtype rows once
    sub_rows = {st: ac.subtype_rows(meta, st) for st in ac.SUBTYPES}

    recs = []
    for ri, R in enumerate(my_regs):
        cols = np.flatnonzero(edge_src == R)
        for st in ac.SUBTYPES:
            rows = sub_rows[st]
            if len(rows) < ac.MIN_CELLS_PER_SECTION:
                continue
            # deterministic per (regulator, subtype) rng -> reproducible permutations
            seed = (ac.SEED, int(np.frombuffer(R.encode()[:4].ljust(4, b"\0"),
                                               dtype=">u4")[0]), ac.SUBTYPES.index(st))
            rng = np.random.default_rng(seed)
            stat = regulator_subtype(X, meta, rows, cols, rng)
            stat.update({"subtype": st, "regulator": R, "d_edges": int(len(cols))})
            recs.append(stat)
        if (ri + 1) % 5 == 0:
            print(f"  ... {ri + 1}/{len(my_regs)} regulators done", flush=True)

    df = pd.DataFrame(recs)
    cols_order = ["subtype", "regulator", "d_edges", "I_R", "E_perm", "p_perm",
                  "I_R_eqsec", "p_eqsec", "n_eqsec",
                  "n_active", "n_cells_used", "n_sections",
                  "I_s_min", "I_s_median", "I_s_max"]
    df = df[cols_order]
    out = SHARD_DIR / f"sq4_regulator_morans_shard{args.shard:03d}.csv"
    df.to_csv(out, index=False)
    print(f"[shard {args.shard}] wrote {len(df)} rows -> {out}", flush=True)


if __name__ == "__main__":
    main()
