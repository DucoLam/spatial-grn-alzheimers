"""
SQ4 / pc_distribution -- where does each top spatial regulator's weight sit across PCs?

"Undo the PCA": each SQ2/SQ3 principal component is a weighted recipe over all
246,890 edges (the SVD loadings). For the top-20 regulators by spatial autocorrelation
(I_R, from the Stage-1 screen), sum each regulator's edge-loadings within every PC to
get how its weight is DISTRIBUTED across the 20 PCs -- and flag the AD-correlated PCs
that SQ3 linked to disease severity.

Key facts (verified):
  * SQ2 and SQ3 both run TruncatedSVD(20, random_state=42) on X[subtype_rows]
    (the same per-subtype edge-weight matrix), deterministic in the fixed container
    -> identical components_ (same values, signs, ordering). So PC indices here match
    SQ3's PC2/3/7 etc. We re-verify by recomputing each PC's Moran's I and matching the
    stored SQ2 table (numbering guard) before trusting any PC label.
  * Loadings are PER EDGE; a regulator's contribution = sum of |loading| over its edges.
  * Per-regulator row-normalisation (each regulator's distribution sums to 1) cancels the
    regulon-size confound for the SHAPE of the distribution (same edges in every PC).

AD-correlated PC leads from SQ3 (1-based): L2/3 IT {2,3,7} . Oligodendrocyte {1} .
Astrocyte {1,3} . L4 IT {8}.

Outputs -> analysis/SQ4/
  sq4_pc_distribution.csv          long: subtype, regulator, I_R, d_edges, pc, frac, is_ad_pc
  sq4_pc_distribution_<slug>.png   heatmap (top-20 regulators x 20 PCs), AD-PC cols boxed
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from sklearn.decomposition import TruncatedSVD

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ4")
SQ2_MORANS = ac.ANALYSIS / "SQ2" / "morans_i"
MASTER = OUT / "sq4_regulator_morans.csv"
TOPN = 20

# SQ3 disease-correlated PC leads (1-based PC numbers)
AD_PCS = {"L2/3 IT": {2, 3, 7}, "Oligodendrocyte": {1},
          "Astrocyte": {1, 3}, "L4 IT": {8}}


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
    edge_src = idx["edge_src"].astype(str)
    if edge_src.shape[0] != X.shape[1]:
        raise ValueError(f"edge_src {edge_src.shape[0]} != n_cols {X.shape[1]}")

    master = pd.read_csv(MASTER)
    recs = []

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < ac.MIN_CELLS_PER_SECTION:
            print(f"[{st}] too few cells, skipping"); continue
        n_comp = min(ac.N_PCS, len(rows) - 1)

        svd = TruncatedSVD(n_components=n_comp, random_state=ac.SEED)
        G = svd.fit_transform(X[rows])             # cell scores (for the numbering guard)
        comps = np.abs(svd.components_)            # (n_comp, n_edges) |loading| per edge

        # --- numbering guard: recompute per-PC Moran's I, match stored SQ2 table ---
        W = ac.within_section_knn_W(meta, rows)
        sq2_csv = SQ2_MORANS / f"morans_{ac.SLUG[st]}.csv"
        if sq2_csv.exists():
            sq2 = pd.read_csv(sq2_csv)
            sq2 = sq2[sq2.space == "GRN"].sort_values("pc")
            ok = True
            for c in range(min(n_comp, len(sq2))):
                I_here = ac.morans_I(G[:, c] - G[:, c].mean(), W)
                I_sq2 = float(sq2.iloc[c]["morans_I"])
                if not np.isfinite(I_here) or abs(I_here - I_sq2) > 1e-3:
                    ok = False
                    print(f"  [{st}] WARN PC{c+1} Moran mismatch: here={I_here:.4f} "
                          f"sq2={I_sq2:.4f}")
            print(f"[{st}] SVD numbering guard vs SQ2: "
                  f"{'OK' if ok else 'MISMATCH -- PC labels suspect'}")
        else:
            print(f"[{st}] no SQ2 morans CSV found ({sq2_csv}); skipping numbering guard")

        # --- top-20 regulators by I_R in this subtype ---
        m = master[(master.subtype == st) & master.I_R.notna()]
        top = m.sort_values("I_R", ascending=False).head(TOPN)["regulator"].tolist()
        ir_map = dict(zip(m.regulator, m.I_R))

        # --- per-regulator loading distribution across PCs ---
        dist = np.zeros((len(top), n_comp))
        d_edges = []
        for i, R in enumerate(top):
            cols = np.flatnonzero(edge_src == R)
            d_edges.append(len(cols))
            c_Rk = comps[:, cols].sum(axis=1)       # (n_comp,) L1 loading of R per PC
            s = c_Rk.sum()
            dist[i] = c_Rk / s if s > 0 else np.nan
            for k in range(n_comp):
                recs.append({"subtype": st, "regulator": R,
                             "I_R": float(ir_map.get(R, np.nan)), "d_edges": len(cols),
                             "pc": k + 1, "frac": float(dist[i, k]),
                             "is_ad_pc": int((k + 1) in AD_PCS.get(st, set()))})

        # --- heatmap (top-PLOTN regulators by I_R) ---
        PLOTN = min(20, len(top))
        distP = dist[:PLOTN]
        topP = top[:PLOTN]
        vmax = np.nanmax(distP) if np.isfinite(np.nanmax(distP)) else 1.0

        fig, ax = plt.subplots(figsize=(12.5, 0.42 * PLOTN + 1.6))
        im = ax.imshow(distP, aspect="auto", cmap="viridis", vmin=0, vmax=vmax)

        # crisp white tile borders
        ax.set_xticks(np.arange(-0.5, n_comp, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, PLOTN, 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=1.4)
        ax.tick_params(which="minor", length=0)

        ax.set_xticks(range(n_comp))
        ax.set_xticklabels([f"PC{k+1}" for k in range(n_comp)], fontsize=11.5)
        ax.set_yticks(range(PLOTN))
        ax.set_yticklabels([f"{R}   $I_R$={ir_map.get(R, np.nan):.3f}" for R in topP],
                           fontsize=11)
        ax.set_xlabel("GRN principal component", fontsize=13.5)

        ax.set_title(f"Loading of top {PLOTN} spatial regulators — {st}",
                     fontsize=15, pad=8)
        for k in range(n_comp):
            if (k + 1) in AD_PCS.get(st, set()):
                ax.add_patch(Rectangle((k - 0.5, -0.5), 1, PLOTN,
                                       fill=False, edgecolor="red", lw=2.6, zorder=5))
                ax.get_xticklabels()[k].set_color("red")
                ax.get_xticklabels()[k].set_fontweight("bold")
        cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
        cbar.set_label("fraction of regulator's total loading", fontsize=11)
        fig.tight_layout()
        fig.savefig(OUT / f"sq4_pc_distribution_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)

        # quick text readout: dominant PC per top regulator + AD-PC weight
        ad = sorted(AD_PCS.get(st, set()))
        ad_idx = [p - 1 for p in ad if p - 1 < n_comp]
        print(f"[{st}] top-{len(top)} regulators | AD PCs {ad}:")
        for i, R in enumerate(top[:10]):
            dom = int(np.nanargmax(dist[i])) + 1
            ad_w = float(dist[i, ad_idx].sum()) if ad_idx else 0.0
            print(f"    {R:<14} I={ir_map.get(R, np.nan):.3f}  "
                  f"dominant=PC{dom}  AD-PC weight={ad_w:.2f}")

    pd.DataFrame(recs).to_csv(OUT / "sq4_pc_distribution.csv", index=False)
    print("pc_distribution done ->", OUT)


if __name__ == "__main__":
    main()
