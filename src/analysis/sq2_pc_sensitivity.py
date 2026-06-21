"""
SQ2 robustness -- sensitivity of Moran's I and clustering purity to the number
of GRN principal components (10 / 20 / 50). 20 is the main-analysis value.

Outputs -> analysis/SQ2/pc_sensitivity/
  pc_sensitivity.csv  subtype, n_pcs, morans_meanI, morans_maxI, morans_sig,
                      cluster_k, purity, purity_null, purity_p
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "pc_sensitivity")
PC_LIST = [10, 20, 50]


def neighbour_purity(labels, W, rng, n_perm=ac.N_PERM):
    Wc = W.tocoo(); a, b = Wc.row, Wc.col
    obs = float(np.mean(labels[a] == labels[b]))
    n = len(labels); perm = np.empty(n_perm)
    for p in range(n_perm):
        perm[p] = np.mean(labels[rng.permutation(n)][a] == labels[rng.permutation(n)][b])
    return obs, float(perm.mean()), (np.sum(perm >= obs) + 1) / (n_perm + 1)


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    out = []
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) < ac.MIN_CELLS_PER_SECTION:
            continue
        W = ac.within_section_knn_W(meta, rows, k=ac.KNN_K)
        for npc in PC_LIST:
            ncomp = min(npc, len(rows) - 1)
            G = TruncatedSVD(n_components=ncomp, random_state=ac.SEED).fit_transform(X[rows])
            Is, nsig = [], 0
            for c in range(ncomp):
                rng = np.random.default_rng(ac.SEED + npc * 1000 + c)
                I, _, p = ac.morans_I_perm(G[:, c], W, rng)
                Is.append(I); nsig += int(p < 0.05)
            # clustering at this PC count
            best = None
            for k in range(2, 9):
                km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
                sil = silhouette_score(G, km.labels_, sample_size=min(len(rows), 2000),
                                       random_state=ac.SEED)
                if best is None or sil > best[1]:
                    best = (k, sil, km.labels_)
            ck, _, labels = best
            rng = np.random.default_rng(ac.SEED + npc)
            pur, null, pp = neighbour_purity(labels, W, rng)
            out.append({"subtype": st, "n_pcs": ncomp,
                        "morans_meanI": float(np.mean(Is)), "morans_maxI": float(np.max(Is)),
                        "morans_sig": nsig, "cluster_k": ck,
                        "purity": pur, "purity_null": null, "purity_p": pp})
            print(f"[{st}] n_pcs={ncomp}: MoranMeanI={np.mean(Is):.3f} sig={nsig}/{ncomp} "
                  f"| k={ck} purity={pur:.3f} null={null:.3f} p={pp:.4f}")
    pd.DataFrame(out).to_csv(OUT / "pc_sensitivity.csv", index=False)
    print("pc_sensitivity done ->", OUT)


if __name__ == "__main__":
    main()
