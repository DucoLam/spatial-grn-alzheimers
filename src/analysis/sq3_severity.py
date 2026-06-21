"""
SQ3 -- relate spatial GRN structure to AD severity (donor-level, exploratory).

For each subtype and each PC count (10/20/50):
  * per donor (>= MIN_DONOR_CELLS cells of that subtype), compute a spatial-
    structure metric from that donor's cells:
       - mean Moran's I across the GRN PCs (within-donor, within-section kNN)
       - neighbour purity (donor's GRN clusters in space)
  * correlate the per-donor metric against the donor's severity (Spearman),
    separately for Braak, Thal, CERAD, ADNC, CPS.

Severity is donor-level (N = donors, ~6-9 per subtype) -> exploratory, underpowered.

Outputs -> analysis/SQ3/
  sq3_perdonor.csv      subtype, n_pcs, donor, n_cells, meanMoranI, purity, + severities
  sq3_correlations.csv  subtype, n_pcs, severity, metric, spearman_r, p, N
  sq3_scatter_<axis>.png (n_pcs=20) per-subtype scatter of meanMoranI vs severity
"""
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from scipy.stats import spearmanr

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ3")
PC_LIST = [10, 20, 50]
MIN_DONOR_CELLS = 100

CERAD = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
ADNC  = {"Not AD": 0, "Low": 1, "Intermediate": 2, "High": 3}

def roman_braak(v):
    m = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}
    return m.get(str(v).replace("Braak", "").strip(), np.nan)

def thal_num(v):
    try:
        return float(str(v).replace("Thal", "").strip())
    except Exception:
        return np.nan


def donor_severity():
    rna = ad.read_h5ad("/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad",
                       backed="r")
    obs = rna.obs.copy()
    obs["Donor ID"] = obs["Donor ID"].astype(str)
    g = obs.groupby("Donor ID", observed=True).first()
    sev = pd.DataFrame(index=g.index)
    sev["CPS"]   = pd.to_numeric(g["Continuous Pseudo-progression Score"], errors="coerce")
    sev["Braak"] = g["Braak"].map(roman_braak)
    sev["Thal"]  = g["Thal"].map(thal_num)
    sev["CERAD"] = g["CERAD score"].astype(str).map(CERAD)
    sev["ADNC"]  = g["Overall AD neuropathological Change"].astype(str).map(ADNC)
    return sev

AXES = ["CPS", "Braak", "Thal", "CERAD", "ADNC"]


def neighbour_purity_val(labels, W):
    Wc = W.tocoo()
    return float(np.mean(labels[Wc.row] == labels[Wc.col]))


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    meta["Donor ID"] = meta["Donor ID"].astype(str)
    sev = donor_severity()

    perdonor, corrs = [], []
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        donors_all = meta["Donor ID"].values[rows]
        for npc in PC_LIST:
            ncomp = min(npc, len(rows) - 1)
            G = TruncatedSVD(n_components=ncomp, random_state=ac.SEED).fit_transform(X[rows])
            # global clusters for purity at this PC count
            best = None
            for k in range(2, 9):
                km = KMeans(n_clusters=k, n_init=10, random_state=ac.SEED).fit(G)
                sil = silhouette_score(G, km.labels_, sample_size=min(len(rows), 2000),
                                       random_state=ac.SEED)
                if best is None or sil > best[1]:
                    best = (k, sil, km.labels_)
            labels = best[2]

            recs = []
            for d in pd.unique(donors_all):
                pos = np.flatnonzero(donors_all == d)        # positions within rows
                if len(pos) < MIN_DONOR_CELLS:
                    continue
                drows = rows[pos]
                Wd = ac.within_section_knn_W(meta, drows, k=ac.KNN_K)
                Gd = G[pos]
                Is = [ac.morans_I(Gd[:, c] - Gd[:, c].mean(), Wd) for c in range(ncomp)]
                meanI = float(np.nanmean(Is))
                pur = neighbour_purity_val(labels[pos], Wd)
                rec = {"subtype": st, "n_pcs": ncomp, "donor": d, "n_cells": int(len(pos)),
                       "meanMoranI": meanI, "purity": pur}
                for ax in AXES:
                    rec[ax] = sev[ax].get(d, np.nan)
                recs.append(rec); perdonor.append(rec)

            df = pd.DataFrame(recs)
            for ax in AXES:
                sub = df.dropna(subset=[ax])
                for metric in ("meanMoranI", "purity"):
                    if len(sub) >= 4 and sub[ax].nunique() > 1:
                        r, p = spearmanr(sub[metric], sub[ax])
                    else:
                        r, p = np.nan, np.nan
                    corrs.append({"subtype": st, "n_pcs": ncomp, "severity": ax,
                                  "metric": metric, "spearman_r": r, "p": p, "N": len(sub)})
            msg = [f"{c['severity']}:r={c['spearman_r']:.2f}(N{c['N']})"
                   for c in corrs if c["subtype"] == st and c["n_pcs"] == ncomp
                   and c["metric"] == "meanMoranI"]
            print(f"[{st}] n_pcs={ncomp}: {len(df)} donors | " + " ".join(msg))

    pd.DataFrame(perdonor).to_csv(OUT / "sq3_perdonor.csv", index=False)
    pd.DataFrame(corrs).to_csv(OUT / "sq3_correlations.csv", index=False)

    # scatter plots at n_pcs=20 (meanMoranI vs each axis)
    pd_df = pd.DataFrame(perdonor)
    d20 = pd_df[pd_df.n_pcs == 20]
    for ax_name in AXES:
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        for axp, st in zip(axes, ac.SUBTYPES):
            s = d20[d20.subtype == st].dropna(subset=[ax_name])
            axp.scatter(s[ax_name], s["meanMoranI"])
            for _, r in s.iterrows():
                axp.annotate(r["donor"].split(".")[-1], (r[ax_name], r["meanMoranI"]), fontsize=7)
            axp.set_title(st); axp.set_xlabel(ax_name); axp.set_ylabel("mean Moran's I")
        fig.suptitle(f"SQ3: per-donor spatial GRN structure vs {ax_name} (n_pcs=20)")
        fig.tight_layout()
        fig.savefig(OUT / f"sq3_scatter_{ax_name}.png", dpi=150)
        plt.close(fig)
    print("sq3 done ->", OUT)


if __name__ == "__main__":
    main()
