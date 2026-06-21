"""
SQ2 robustness -- sensitivity to the spatial-neighbour k (the k=6 choice).

Re-runs Moran's I (mean/max I, # significant PCs) and neighbour purity for
spatial kNN k in {4,6,8,10,15}, per subtype. Cluster labels are fixed once per
subtype (silhouette-optimal k-means on GRN PCs) and reused across all spatial k,
so ONLY the spatial-neighbour size varies. If the result holds across k, the
k=6 choice in the main analysis is not driving it.

Outputs -> analysis/SQ2/spatialk_sweep/
  spatialk_sweep.csv   subtype, spatial_k, morans meanI/maxI/#sig, purity/null/p
  spatialk_sweep.png   mean Moran's I and purity vs spatial k, per subtype
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

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "spatialk_sweep")
SPATIAL_KS = [4, 6, 8, 10, 15]
CLUSTER_KRANGE = range(2, 9)


def neighbour_purity(labels, W, rng, n_perm=ac.N_PERM):
    Wc = W.tocoo(); a, b = Wc.row, Wc.col
    obs = float(np.mean(labels[a] == labels[b]))
    n = len(labels); perm = np.empty(n_perm)
    for p in range(n_perm):
        lp = labels[rng.permutation(n)]
        perm[p] = np.mean(lp[a] == lp[b])
    return obs, float(perm.mean()), (np.sum(perm >= obs) + 1) / (n_perm + 1)


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    out = []
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < ac.MIN_CELLS_PER_SECTION:
            continue
        n_comp = min(ac.N_PCS, len(rows) - 1)
        G = TruncatedSVD(n_components=n_comp, random_state=ac.SEED).fit_transform(X[rows])

        # fix cluster labels once (silhouette-optimal), reuse across spatial k
        best = None
        for k in CLUSTER_KRANGE:
            km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
            sil = silhouette_score(G, km.labels_, sample_size=min(len(rows), 2000),
                                   random_state=ac.SEED)
            if best is None or sil > best[1]:
                best = (k, sil, km.labels_)
        ck, _, labels = best

        for sk in SPATIAL_KS:
            W = ac.within_section_knn_W(meta, rows, k=sk)
            Is, nsig = [], 0
            for c in range(n_comp):
                rng = np.random.default_rng(ac.SEED + sk * 100 + c)
                I, _, p = ac.morans_I_perm(G[:, c], W, rng)
                Is.append(I); nsig += int(p < 0.05)
            rng = np.random.default_rng(ac.SEED + sk)
            pur, null, pp = neighbour_purity(labels, W, rng)
            out.append({"subtype": st, "spatial_k": sk, "cluster_k": ck,
                        "morans_meanI": float(np.mean(Is)), "morans_maxI": float(np.max(Is)),
                        "morans_sig": nsig, "n_pcs": n_comp,
                        "purity": pur, "purity_null": null, "purity_p": pp})
            print(f"[{st}] spatial_k={sk}: MoranMeanI={np.mean(Is):.3f} "
                  f"sig={nsig}/{n_comp} | purity={pur:.3f} null={null:.3f} p={pp:.4f}")

    df = pd.DataFrame(out)
    df.to_csv(OUT / "spatialk_sweep.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for st in ac.SUBTYPES:
        d = df[df.subtype == st]
        if len(d) == 0:
            continue
        axes[0].plot(d.spatial_k, d.morans_meanI, "-o", label=st)
        axes[1].plot(d.spatial_k, d.purity, "-o", label=st)
        axes[1].plot(d.spatial_k, d.purity_null, "--", color="0.7", lw=1)
    axes[0].set_title("Mean Moran's I vs spatial k"); axes[0].set_xlabel("spatial kNN k")
    axes[0].set_ylabel("mean Moran's I"); axes[0].axvline(6, color="k", lw=0.6, ls=":")
    axes[1].set_title("Neighbour purity vs spatial k (dashed = chance)")
    axes[1].set_xlabel("spatial kNN k"); axes[1].set_ylabel("purity")
    axes[1].axvline(6, color="k", lw=0.6, ls=":")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("SQ2 robustness to spatial-neighbour k (k=6 = main analysis)")
    fig.tight_layout()
    fig.savefig(OUT / "spatialk_sweep.png", dpi=160)
    plt.close(fig)
    print("spatialk_sweep done ->", OUT)


if __name__ == "__main__":
    main()
