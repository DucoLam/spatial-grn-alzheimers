"""
SQ2 / grn_clustering robustness -- neighbour purity across ALL k.

Addresses the "k was cherry-picked" concern: instead of trusting the single
silhouette-optimal k, we run the spatial co-location test (neighbour purity +
permutation p) for every k in 2..8, per subtype, and show the result holds
regardless of k. k is chosen by silhouette (GRN-only, spatially blind); the
spatial test is independent of that choice.

Outputs -> analysis/SQ2/grn_clustering/
  purity_across_k.csv   subtype, k, silhouette, purity, purity_null, purity_p
  purity_across_k.png   purity vs null across k, per subtype (2x2)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "grn_clustering")
K_RANGE = list(range(2, 9))


def neighbour_purity(labels, W, rng, n_perm=ac.N_PERM):
    Wc = W.tocoo(); a, b = Wc.row, Wc.col
    obs = float(np.mean(labels[a] == labels[b]))
    n = len(labels)
    perm = np.empty(n_perm)
    for p in range(n_perm):
        lp = labels[rng.permutation(n)]
        perm[p] = np.mean(lp[a] == lp[b])
    pval = (np.sum(perm >= obs) + 1) / (n_perm + 1)
    return obs, float(perm.mean()), pval


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    rows_out = []

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < max(ac.MIN_CELLS_PER_SECTION, 2 * max(K_RANGE)):
            print(f"[{st}] too few cells, skipping"); continue
        n_comp = min(ac.N_PCS, len(rows) - 1)
        G = TruncatedSVD(n_components=n_comp, random_state=ac.SEED).fit_transform(X[rows])
        W = ac.within_section_knn_W(meta, rows)
        sil_n = min(len(rows), 2000)

        for k in K_RANGE:
            km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
            sil = silhouette_score(G, km.labels_, sample_size=sil_n, random_state=ac.SEED)
            rng = np.random.default_rng(ac.SEED + k)
            pur, null, pval = neighbour_purity(km.labels_, W, rng)
            rows_out.append({"subtype": st, "k": k, "silhouette": float(sil),
                             "purity": pur, "purity_null": null, "purity_p": pval})
            print(f"[{st}] k={k}: sil={sil:.3f} purity={pur:.3f} "
                  f"null={null:.3f} p={pval:.4f}")

    df = pd.DataFrame(rows_out)
    df.to_csv(OUT / "purity_across_k.csv", index=False)

    # figure: purity vs null across k, per subtype
    subs = [s for s in ac.SUBTYPES if s in df.subtype.values]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
    for ax, st in zip(axes.ravel(), subs):
        d = df[df.subtype == st].sort_values("k")
        ax.plot(d.k, d.purity, "-o", color="#1f77b4", label="observed purity")
        ax.plot(d.k, d.purity_null, "--s", color="#999999", label="chance (null)")
        for _, r in d.iterrows():
            if r.purity_p < 0.05:
                ax.text(r.k, r.purity + 0.01, "*", ha="center", color="#1f77b4")
        ax.set_title(st); ax.set_ylabel("neighbour purity"); ax.set_xlabel("k")
    axes.ravel()[0].legend(frameon=False, fontsize=8)
    fig.suptitle("SQ2: spatial neighbour purity vs chance, across all k "
                 "(* = p<0.05, permutation)")
    fig.tight_layout()
    fig.savefig(OUT / "purity_across_k.png", dpi=160)
    plt.close(fig)
    print("purity_across_k done ->", OUT)


if __name__ == "__main__":
    main()
