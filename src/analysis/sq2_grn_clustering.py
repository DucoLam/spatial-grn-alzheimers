"""
SQ2 / grn_clustering -- connect the result back to ScReNI's own method.

Per subtype: cluster cells purely by their GRN (KMeans on GRN PCs; k chosen by
silhouette over 2..8), with NO spatial information used. Then ask whether those
GRN-defined clusters are spatially organised:
  * neighbour purity -- fraction of within-section kNN pairs that share a GRN
    cluster, vs a label-permutation null (spatial co-location of GRN clusters).
  * per-section Kruskal-Wallis of x and y by cluster (descriptive: do clusters
    sit at different positions within a section). No depth-from-pia axis exists,
    so x and y are reported separately rather than a single depth call.

Outputs -> analysis/SQ2/grn_clustering/
  clustering_summary.csv          per subtype: k, silhouette, purity, purity_p
  clustering_labels_<slug>.csv    cell_id, section, x, y, grn_cluster
  clustering_<slug>.png           per-section scatter coloured by GRN cluster
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import kruskal

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "grn_clustering")
K_RANGE = range(2, 9)


def neighbour_purity(labels, W, rng, n_perm=ac.N_PERM):
    """Fraction of kNN edges whose endpoints share a cluster, vs permuted null."""
    Wc = W.tocoo()
    a, b = Wc.row, Wc.col
    obs = float(np.mean(labels[a] == labels[b]))
    perm = np.empty(n_perm)
    n = len(labels)
    for p in range(n_perm):
        lp = labels[rng.permutation(n)]
        perm[p] = np.mean(lp[a] == lp[b])
    p_val = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    return obs, float(perm.mean()), p_val


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    summary = []

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < max(ac.MIN_CELLS_PER_SECTION, 2 * max(K_RANGE)):
            print(f"[{st}] too few cells, skipping"); continue

        n_comp = min(ac.N_PCS, len(rows) - 1)
        G = TruncatedSVD(n_components=n_comp, random_state=ac.SEED).fit_transform(X[rows])

        # pick k by silhouette (sampled for speed on large subtypes)
        sil_n = min(len(rows), 2000)
        best = None
        for k in K_RANGE:
            km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
            sil = silhouette_score(G, km.labels_, sample_size=sil_n,
                                   random_state=ac.SEED)
            if best is None or sil > best[1]:
                best = (k, sil, km.labels_)
        k, sil, labels = best

        # spatial organisation of the GRN clusters
        W = ac.within_section_knn_W(meta, rows)
        rng = np.random.default_rng(ac.SEED)
        purity, purity_null, purity_p = neighbour_purity(labels, W, rng)

        # per-section Kruskal-Wallis on x and y (descriptive)
        sect = meta["section"].values[rows]
        xy = meta[["x", "y"]].values[rows]
        kw_x, kw_y, nsec_tested = [], [], 0
        for s in pd.unique(sect):
            m = sect == s
            lab_s = labels[m]
            if len(np.unique(lab_s)) < 2 or m.sum() < ac.MIN_CELLS_PER_SECTION:
                continue
            groups_x = [xy[m, 0][lab_s == u] for u in np.unique(lab_s)]
            groups_y = [xy[m, 1][lab_s == u] for u in np.unique(lab_s)]
            try:
                kw_x.append(kruskal(*groups_x).pvalue)
                kw_y.append(kruskal(*groups_y).pvalue)
                nsec_tested += 1
            except ValueError:
                pass
        frac_sig = np.nan
        if nsec_tested:
            sig = np.mean((np.array(kw_x) < 0.05) | (np.array(kw_y) < 0.05))
            frac_sig = float(sig)

        summary.append({"subtype": st, "n_cells": int(len(rows)), "k": k,
                        "silhouette": float(sil),
                        "neighbour_purity": purity, "purity_null": purity_null,
                        "purity_p": purity_p,
                        "sections_tested": nsec_tested,
                        "frac_sections_kw_sig": frac_sig})
        print(f"[{st}] k={k} sil={sil:.3f} | purity {purity:.3f} "
              f"(null {purity_null:.3f}, p={purity_p:.4f}) | "
              f"KW-sig sections {frac_sig}")

        # save labels
        lab_df = pd.DataFrame({"cell_id": cell_ids[rows], "section": sect,
                               "x": xy[:, 0], "y": xy[:, 1], "grn_cluster": labels})
        lab_df.to_csv(OUT / f"clustering_labels_{ac.SLUG[st]}.csv", index=False)

        # per-section scatter (up to 6 largest sections)
        top = lab_df["section"].value_counts().head(6).index.tolist()
        ncol = min(3, len(top)); nrow = int(np.ceil(len(top) / ncol)) if top else 1
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.6 * nrow),
                                 squeeze=False)
        for ax in axes.ravel():
            ax.axis("off")
        for i, s in enumerate(top):
            ax = axes.ravel()[i]; ax.axis("on")
            d = lab_df[lab_df.section == s]
            ax.scatter(d.x, d.y, c=d.grn_cluster, cmap="tab10", s=6)
            ax.set_title(f"section {s} (n={len(d)})", fontsize=9)
            ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        fig.suptitle(f"SQ2 GRN clusters in space — {st} (k={k})")
        fig.tight_layout()
        fig.savefig(OUT / f"clustering_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)

    pd.DataFrame(summary).to_csv(OUT / "clustering_summary.csv", index=False)
    print("grn_clustering done ->", OUT)


if __name__ == "__main__":
    main()
