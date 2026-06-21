"""
SQ2 / conditional -- does GRN spatial structure survive controlling for expression?

This is the circularity-control analysis. Tangram places cells by expression
similarity, so expression is spatially autocorrelated by construction and the GRN
(derived from expression) may inherit that. Two conditional tests ask whether GRN
spatial structure exceeds what expression alone explains:

  (1) Partial Mantel:  r(GRN-dissim, spatial-dist | expression-dissim), pooled over
      within-section pairs, permutation p (GRN labels shuffled within section).
  (2) Expression-residualised Moran's I:  regress each GRN principal component on
      ALL expression PCs, then test whether the RESIDUAL still clusters in space.
      If residual Moran's I stays positive & significant, the spatial structure is
      not merely inherited from expression.

Outputs -> analysis/SQ2/conditional/
  partial_mantel.csv          per subtype: r_GS, r_GE, r_SE, partial_r, p_perm
  residual_morans_<slug>.csv  per PC: morans_raw, morans_resid, p_resid, expr_R2
  residual_morans.png         raw vs residual Moran's I, 2x2 subtypes
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.linear_model import LinearRegression

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "conditional")


# --------------------------------------------------------------------------
# (1) partial Mantel
# --------------------------------------------------------------------------
def tri_blocks(X, expr, meta, rows, rng):
    """Per-section spatial / GRN / expression dissimilarities on the SAME cells."""
    coords = meta[["x", "y"]].values
    sect = meta["section"].values
    blocks = []
    for s in pd.unique(sect[rows]):
        sel = rows[sect[rows] == s]
        if len(sel) < ac.MIN_CELLS_PER_SECTION:
            continue
        if len(sel) > ac.MAX_CELLS_PER_SECTION:
            sel = np.sort(rng.choice(sel, ac.MAX_CELLS_PER_SECTION, replace=False))
        xy = coords[sel]
        diff = xy[:, None, :] - xy[None, :, :]
        S = np.sqrt((diff ** 2).sum(-1)).astype(np.float32)
        Gsq = ac._cosine_dissim(X[sel])          # square (for permutation)
        Esq = ac._cosine_dissim(expr[sel])
        iu = np.triu_indices(len(sel), k=1)
        blocks.append({"n": len(sel), "iu": iu,
                       "S": S[iu], "E": Esq[iu], "Gsq": Gsq, "G": Gsq[iu]})
    return blocks


def partial_r(r_gs, r_ge, r_se):
    denom = np.sqrt((1 - r_ge ** 2) * (1 - r_se ** 2))
    return (r_gs - r_ge * r_se) / denom if denom > 0 else np.nan


def partial_mantel(blocks, rng, n_perm=ac.N_PERM):
    S = np.concatenate([b["S"] for b in blocks])
    E = np.concatenate([b["E"] for b in blocks])
    G = np.concatenate([b["G"] for b in blocks])
    r_gs, r_ge, r_se = ac._pearson(G, S), ac._pearson(G, E), ac._pearson(S, E)
    obs = partial_r(r_gs, r_ge, r_se)
    perm = np.empty(n_perm)
    for p in range(n_perm):
        Gp = []
        for b in blocks:
            o = rng.permutation(b["n"])
            Gp.append(b["Gsq"][np.ix_(o, o)][b["iu"]])
        Gp = np.concatenate(Gp)
        perm[p] = partial_r(ac._pearson(Gp, S), ac._pearson(Gp, E), r_se)
    p_val = (np.sum(np.abs(perm) >= abs(obs)) + 1) / (n_perm + 1)
    return {"r_GS": r_gs, "r_GE": r_ge, "r_SE": r_se,
            "partial_r": obs, "p_perm": p_val,
            "n_pairs": int(len(S)), "n_sections": len(blocks)}


# --------------------------------------------------------------------------
# (2) expression-residualised Moran's I
# --------------------------------------------------------------------------
def residual_morans(X, expr, meta, rows):
    W = ac.within_section_knn_W(meta, rows)
    n_comp = min(ac.N_PCS, len(rows) - 1)
    Gpc = TruncatedSVD(n_components=n_comp, random_state=ac.SEED).fit_transform(X[rows])
    Epc = PCA(n_components=min(n_comp, expr.shape[1]),
              random_state=ac.SEED).fit_transform(expr[rows])
    recs = []
    for c in range(n_comp):
        y = Gpc[:, c]
        lr = LinearRegression().fit(Epc, y)
        resid = y - lr.predict(Epc)
        r2 = lr.score(Epc, y)              # how much of this GRN PC expression explains
        rng = np.random.default_rng(ac.SEED + c)
        I_raw, _, p_raw = ac.morans_I_perm(y, W, rng)
        rng = np.random.default_rng(ac.SEED + 5000 + c)
        I_res, _, p_res = ac.morans_I_perm(resid, W, rng)
        recs.append({"pc": c + 1, "expr_R2": float(r2),
                     "morans_raw": I_raw, "p_raw": p_raw,
                     "morans_resid": I_res, "p_resid": p_res})
    return pd.DataFrame(recs)


def main():
    X, cell_ids, meta, expr = ac.load_cache(with_expr=True)
    pm_rows, resid_tables = [], {}

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < ac.MIN_CELLS_PER_SECTION:
            print(f"[{st}] too few cells, skipping"); continue

        # (1) partial Mantel
        rng = np.random.default_rng(ac.SEED)
        blocks = tri_blocks(X, expr, meta, rows, rng)
        if blocks:
            pm = partial_mantel(blocks, rng)
            pm.update({"subtype": st, "n_cells": int(len(rows))})
            pm_rows.append(pm)
            print(f"[{st}] partial Mantel: r(GS)={pm['r_GS']:.4f} "
                  f"r(GS|E)={pm['partial_r']:.4f} p={pm['p_perm']:.4f}  "
                  f"(r_GE={pm['r_GE']:.3f}, r_SE={pm['r_SE']:.3f})")

        # (2) residualised Moran's I
        df = residual_morans(X, expr, meta, rows)
        df.to_csv(OUT / f"residual_morans_{ac.SLUG[st]}.csv", index=False)
        resid_tables[st] = df
        n_raw = int((df.p_raw < 0.05).sum())
        n_res = int((df.p_resid < 0.05).sum())
        print(f"[{st}] Moran's I sig PCs: raw {n_raw}/{len(df)} -> "
              f"residual {n_res}/{len(df)} | "
              f"mean I {df.morans_raw.mean():.3f} -> {df.morans_resid.mean():.3f} "
              f"| mean expr-R2 {df.expr_R2.mean():.2f}")

    if pm_rows:
        cols = ["subtype", "n_cells", "n_sections", "n_pairs",
                "r_GS", "r_GE", "r_SE", "partial_r", "p_perm"]
        pd.DataFrame(pm_rows)[cols].to_csv(OUT / "partial_mantel.csv", index=False)

    # figure: raw vs residual Moran's I, 2x2
    if resid_tables:
        fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
        for ax, (st, df) in zip(axes.ravel(), resid_tables.items()):
            x = df["pc"].values
            ax.bar(x - 0.2, df["morans_raw"], width=0.4, color="#1f77b4",
                   label="raw GRN")
            ax.bar(x + 0.2, df["morans_resid"], width=0.4, color="#d62728",
                   label="residual (expression removed)")
            ax.axhline(0, color="k", lw=0.8)
            ax.set_title(st); ax.set_xticks(x[::2]); ax.set_ylabel("Moran's I")
        axes.ravel()[0].legend(frameon=False)
        for ax in axes[1]:
            ax.set_xlabel("GRN principal component")
        fig.suptitle("SQ2 conditional: GRN spatial autocorrelation before vs after "
                     "removing expression", y=1.0)
        fig.tight_layout()
        fig.savefig(OUT / "residual_morans.png", dpi=160)
        plt.close(fig)
    print("conditional done ->", OUT)


if __name__ == "__main__":
    main()
