"""
SQ3 -- cluster-composition vs AD severity (donor-level, exploratory).

Per subtype: cluster cells into k=3 GRN clusters (KMeans on 20 GRN PCs). For each
donor (>= MIN_DONOR_CELLS cells), compute the fraction of its cells in each cluster.
Then correlate (Spearman, across donors) each cluster's per-donor proportion against
each severity axis:
  Braak, Thal, CERAD, ADNC, CPS, severity-PC1, severity-PC2, tau-minus-amyloid.
Benjamini-Hochberg FDR across the whole grid. N = donors (~6-9) -> exploratory.

Outputs -> analysis/SQ3/
  sq3_composition_corr.csv   subtype, cluster, severity, spearman_r, p, fdr_q, N
  sq3_composition_props.csv  per-donor cluster proportions + severities
  sq3_severity_loadings.csv  severity-PCA loadings (what PC1/PC2 mean)
  sq3_composition_<slug>.png heatmap clusters x severity (r, * = FDR<0.05)
"""
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD, PCA
from sklearn.cluster import KMeans
from scipy.stats import spearmanr

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ3")
# per-subtype k from the silhouette peaks (L2/3 clean at 2, L4 at 3, Astro local peak at 4;
# Oligo has no real peak -> k=3 as an interpretable default)
K_PER_SUBTYPE = {"Astrocyte": 4, "L2/3 IT": 2, "L4 IT": 3, "Oligodendrocyte": 3}
N_PCS = 20
MIN_DONOR_CELLS = 100

CERAD = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
ADNC  = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}
def braak_n(v): return {"0":0,"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6}.get(str(v).replace("Braak","").strip(), np.nan)
def thal_n(v):
    try: return float(str(v).replace("Thal","").strip())
    except Exception: return np.nan

PATH_AXES = ["Braak", "Thal", "CERAD", "ADNC", "CPS"]


def bh_fdr(pvals):
    p = np.asarray(pvals, float); ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    pv = p[ok]; n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = np.empty(n); ranked[order] = np.arange(1, n + 1)
    qv = pv * n / ranked
    qs = qv[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qs, 0, 1)
    q[ok] = out
    return q


def donor_severity():
    rna = ad.read_h5ad("/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad", backed="r")
    obs = rna.obs.copy(); obs["Donor ID"] = obs["Donor ID"].astype(str)
    g = obs.groupby("Donor ID", observed=True).first()
    s = pd.DataFrame(index=g.index)
    s["CPS"]   = pd.to_numeric(g["Continuous Pseudo-progression Score"], errors="coerce")
    s["Braak"] = g["Braak"].map(braak_n)
    s["Thal"]  = g["Thal"].map(thal_n)
    s["CERAD"] = g["CERAD score"].astype(str).map(CERAD)
    s["ADNC"]  = g["Overall AD neuropathological Change"].astype(str).map(ADNC)
    s = s.apply(lambda col: pd.to_numeric(col, errors="coerce"))  # kill categorical dtype
    return s


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    meta["Donor ID"] = meta["Donor ID"].astype(str)
    sev = donor_severity()

    # severity structure: standardize complete donors, PCA + tau-amyloid contrast
    sev_c = sev.dropna()
    Z = (sev_c - sev_c.mean()) / sev_c.std(ddof=0)
    pca = PCA(n_components=2, random_state=ac.SEED).fit(Z.values)
    pcs = pca.transform(Z.values)
    sev2 = sev.copy()
    sev2["sevPC1"] = np.nan; sev2["sevPC2"] = np.nan
    sev2.loc[sev_c.index, "sevPC1"] = pcs[:, 0]
    sev2.loc[sev_c.index, "sevPC2"] = pcs[:, 1]
    sev2["tau_minus_amyloid"] = (Z["Braak"] - Z["Thal"]).reindex(sev2.index)
    load = pd.DataFrame(pca.components_.T, index=sev_c.columns, columns=["sevPC1", "sevPC2"])
    load.to_csv(OUT / "sq3_severity_loadings.csv")
    print("severity-PCA loadings:\n", load[["sevPC1", "sevPC2"]].round(3))
    print("var explained PC1,PC2:", np.round(pca.explained_variance_ratio_, 3))
    AXES = PATH_AXES + ["sevPC1", "sevPC2", "tau_minus_amyloid"]

    props_rows, corr_rows = [], []
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        k_st = K_PER_SUBTYPE[st]
        ncomp = min(N_PCS, len(rows) - 1)
        G = TruncatedSVD(n_components=ncomp, random_state=ac.SEED).fit_transform(X[rows])
        labels = KMeans(n_clusters=k_st, n_init=10, random_state=ac.SEED).fit_predict(G)
        donors = meta["Donor ID"].values[rows]
        # per-donor cluster proportions
        recs = {}
        for d in pd.unique(donors):
            m = donors == d
            if m.sum() < MIN_DONOR_CELLS:
                continue
            lab = labels[m]
            props = [(lab == c).mean() for c in range(k_st)]
            row = {"subtype": st, "donor": d, "n_cells": int(m.sum())}
            for c in range(k_st):
                row[f"prop_c{c}"] = props[c]
            for ax in AXES:
                row[ax] = sev2[ax].get(d, np.nan)
            recs[d] = row; props_rows.append(row)
        df = pd.DataFrame(recs).T
        for c in range(k_st):
            for ax in AXES:
                sub = df.dropna(subset=[ax, f"prop_c{c}"])
                if len(sub) >= 4 and sub[ax].astype(float).nunique() > 1 and sub[f"prop_c{c}"].astype(float).nunique() > 1:
                    r, p = spearmanr(sub[f"prop_c{c}"].astype(float), sub[ax].astype(float))
                else:
                    r, p = np.nan, np.nan
                corr_rows.append({"subtype": st, "cluster": f"c{c}", "severity": ax,
                                  "spearman_r": r, "p": p, "N": len(sub)})
        print(f"[{st}] k={k_st}, {df.shape[0]} donors with >= {MIN_DONOR_CELLS} cells")

    cdf = pd.DataFrame(corr_rows)
    cdf["fdr_q"] = bh_fdr(cdf["p"].values)
    cdf.to_csv(OUT / "sq3_composition_corr.csv", index=False)
    pd.DataFrame(props_rows).to_csv(OUT / "sq3_composition_props.csv", index=False)

    # heatmaps per subtype
    for st in ac.SUBTYPES:
        d = cdf[cdf.subtype == st]
        if d["spearman_r"].notna().sum() == 0:
            continue
        mat = d.pivot(index="cluster", columns="severity", values="spearman_r").reindex(columns=AXES)
        qmat = d.pivot(index="cluster", columns="severity", values="fdr_q").reindex(columns=AXES)
        fig, ax = plt.subplots(figsize=(9, 3))
        im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(AXES))); ax.set_xticklabels(AXES, rotation=45, ha="right")
        ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                r = mat.values[i, j]; q = qmat.values[i, j]
                if not np.isnan(r):
                    star = "*" if (not np.isnan(q) and q < 0.05) else ""
                    ax.text(j, i, f"{r:.2f}{star}", ha="center", va="center", fontsize=8)
        ax.set_title(f"SQ3 cluster-composition vs severity — {st} (k={K_PER_SUBTYPE[st]}, *=FDR<0.05)")
        fig.colorbar(im, label="Spearman r")
        fig.tight_layout()
        fig.savefig(OUT / f"sq3_composition_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)
    print("sq3_composition done ->", OUT)


if __name__ == "__main__":
    main()
