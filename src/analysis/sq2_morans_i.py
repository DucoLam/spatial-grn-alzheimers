"""
SQ2 / morans_i -- which GRN features are spatially structured?

Per subtype: reduce the GRN edge-weight matrix to N_PCS principal components
(TruncatedSVD on the sparse matrix), then compute global Moran's I for each PC
using a within-section kNN spatial-weights matrix (neighbours never cross a
section frame). Moran's I > 0 with small p => that GRN component is spatially
autocorrelated (clusters in space rather than being randomly arranged).

The same is done on expression PCs as the circularity baseline.

Outputs -> analysis/SQ2/morans_i/
  morans_<slug>.csv     per-PC: space, pc, var_explained, morans_I, E_perm, p_perm
  morans_summary.png    bar chart of Moran's I per PC (GRN), starred if p<0.05
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD, PCA

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "morans_i")


def main():
    X, cell_ids, meta, expr = ac.load_cache(with_expr=True)

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < ac.MIN_CELLS_PER_SECTION:
            print(f"[{st}] too few cells, skipping"); continue

        W = ac.within_section_knn_W(meta, rows)
        n_comp = min(ac.N_PCS, len(rows) - 1)
        recs = []

        # GRN PCs (sparse -> TruncatedSVD)
        svd = TruncatedSVD(n_components=n_comp, random_state=ac.SEED)
        G = svd.fit_transform(X[rows])
        for c in range(n_comp):
            rng = np.random.default_rng(ac.SEED + c)
            I, E, p = ac.morans_I_perm(G[:, c], W, rng)
            recs.append({"space": "GRN", "pc": c + 1,
                         "var_explained": float(svd.explained_variance_ratio_[c]),
                         "morans_I": I, "E_perm": E, "p_perm": p})

        # expression PCs (dense)
        pca = PCA(n_components=min(n_comp, expr.shape[1]), random_state=ac.SEED)
        Epc = pca.fit_transform(expr[rows])
        for c in range(Epc.shape[1]):
            rng = np.random.default_rng(ac.SEED + 1000 + c)
            I, E, p = ac.morans_I_perm(Epc[:, c], W, rng)
            recs.append({"space": "expression", "pc": c + 1,
                         "var_explained": float(pca.explained_variance_ratio_[c]),
                         "morans_I": I, "E_perm": E, "p_perm": p})

        df = pd.DataFrame(recs)
        df.to_csv(OUT / f"morans_{ac.SLUG[st]}.csv", index=False)
        g = df[df.space == "GRN"]
        nsig = int((g.p_perm < 0.05).sum())
        print(f"[{st}] {len(rows):,} cells | GRN PCs spatially sig (p<0.05): "
              f"{nsig}/{len(g)} | max I = {g.morans_I.max():.3f}")

        # bar chart of GRN Moran's I per PC
        fig, ax = plt.subplots(figsize=(7, 4))
        colors = ["C0" if pp < 0.05 else "0.7" for pp in g.p_perm]
        ax.bar(g.pc, g.morans_I, color=colors)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("GRN principal component")
        ax.set_ylabel("Moran's I")
        ax.set_title(f"SQ2 Moran's I — {st} (filled = p<0.05)")
        fig.tight_layout()
        fig.savefig(OUT / f"morans_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)

    print("morans_i done ->", OUT)


if __name__ == "__main__":
    main()
