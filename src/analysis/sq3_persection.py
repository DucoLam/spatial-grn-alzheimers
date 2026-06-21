"""
SQ3 per-SECTION validation of the per-PC per-donor result.

For each subtype, reuse the SQ2 top-5 spatial PCs (severity-blind). Compute
Moran's I per PC at TWO levels:
  - per DONOR  (reloaded from sq3_perpc_morans.csv -> the validated result)
  - per SECTION (computed here; each MERFISH section = its own coord frame)

Then, per top-5 PC:
  - donor-level Spearman(Moran's I, CPS)        [parametric p, as before]
  - section-level Spearman(Moran's I, CPS)       [observed r]
      significance via DONOR-LEVEL permutation (shuffle donor->CPS, propagate
      to that donor's sections) -> respects pseudo-replication (N_donors is the
      true unit), 9999 perms. BH-FDR within the top-5 set.

Overlay scatter per subtype: faded section points + solid donor points, one
colour per donor, with donor and section trendlines.

Outputs -> analysis/SQ3/
  sq3_persection_morans.csv   subtype, donor, section, pc, morans_I, CPS
  sq3_persection_corr.csv     subtype, pc, pooled_I, donor_r, donor_p,
                              sec_r, sec_perm_p, sec_fdr_q, n_sections, n_donors
  sq3_overlay_<slug>.png      donor+section overlay scatter (top-5 PCs)
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import scipy.sparse as sp
from scipy.stats import spearmanr

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ3")
N_PCS = ac.N_PCS
KNN_K = ac.KNN_K
TOP_K = 5
MIN_SECTION_CELLS = 50      # min cells of the subtype in a section to use it
N_PERM = 9999
SEED = ac.SEED

CERAD_MAP = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}


def donor_cps():
    import anndata as ad
    rna = ad.read_h5ad(
        "/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad",
        backed="r")
    obs = rna.obs.copy(); obs["Donor ID"] = obs["Donor ID"].astype(str)
    g = obs.groupby("Donor ID", observed=True).first()
    return pd.to_numeric(g["Continuous Pseudo-progression Score"], errors="coerce")


def knn_W(coords, sections, k=KNN_K):
    """Row-standardised within-section 6-NN Euclidean weights for a cell subset."""
    n = len(coords)
    rptr, cidx = [], []
    for s in pd.unique(sections):
        loc = np.flatnonzero(sections == s)
        if len(loc) <= 1:
            continue
        kk = min(k, len(loc) - 1)
        nn = NearestNeighbors(n_neighbors=kk + 1).fit(coords[loc])
        _, nbr = nn.kneighbors(coords[loc])
        for a in range(len(loc)):
            for b in nbr[a, 1:]:
                rptr.append(loc[a]); cidx.append(loc[b])
    if not rptr:
        return None
    W = sp.csr_matrix((np.ones(len(rptr)), (rptr, cidx)), shape=(n, n))
    rs = np.asarray(W.sum(1)).ravel(); rs[rs == 0] = 1.0
    return sp.diags(1.0 / rs) @ W


def morans_vec(G, W):
    """Moran's I for each column of G using weights W -> array (n_pcs,)."""
    out = np.full(G.shape[1], np.nan)
    if W is None:
        return out
    Wsum = W.sum(); n = G.shape[0]
    if Wsum == 0:
        return out
    for c in range(G.shape[1]):
        z = G[:, c] - G[:, c].mean()
        den = float(z @ z)
        if den == 0:
            continue
        out[c] = (n / Wsum) * float(z @ (W @ z)) / den
    return out


def bh_fdr(pvals):
    p = np.asarray(pvals, float); ok = ~np.isnan(p)
    q = np.full(len(p), np.nan); pv = p[ok]; n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv); ranked = np.empty(n); ranked[order] = np.arange(1, n + 1)
    qs = (pv * n / ranked)[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qs, 0, 1); q[ok] = out
    return q


def donor_perm_p(sec_I, sec_donor, donor_cps_map, obs_r, rng, nperm=N_PERM):
    """Permute CPS at the DONOR level, propagate to sections, null Spearman."""
    donors = list(donor_cps_map.keys())
    cps_vals = np.array([donor_cps_map[d] for d in donors], float)
    # map each section to its donor index
    didx = np.array([donors.index(d) for d in sec_donor])
    count = 0
    for _ in range(nperm):
        perm = rng.permutation(cps_vals)
        sec_perm_cps = perm[didx]
        r, _ = spearmanr(sec_I, sec_perm_cps)
        if not np.isnan(r) and abs(r) >= abs(obs_r):
            count += 1
    return (count + 1) / (nperm + 1)


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    meta = meta.copy(); meta["Donor ID"] = meta["Donor ID"].astype(str)
    cps = donor_cps()

    # reuse the validated per-donor result for overlay + donor-level corr
    perdonor = pd.read_csv(OUT / "sq3_perpc_morans.csv")  # subtype,donor,pc,morans_I,pooled_I,sq2_sig,CPS

    rng = np.random.default_rng(SEED)
    sec_rows = []
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        meta_st = meta.iloc[rows].reset_index(drop=True)
        n_comp = min(N_PCS, len(rows) - 1)
        G_all = TruncatedSVD(n_components=n_comp, random_state=SEED).fit_transform(X[rows])
        sects = meta_st["section"].values
        donors = meta_st["Donor ID"].values
        coords = meta_st[["x", "y"]].values
        for s in pd.unique(sects):
            m = sects == s
            if m.sum() < MIN_SECTION_CELLS:
                continue
            Is = morans_vec(G_all[m], knn_W(coords[m], sects[m]))
            donor = donors[m][0]
            for c in range(len(Is)):
                sec_rows.append({"subtype": st, "donor": donor, "section": s,
                                 "pc": c + 1, "morans_I": float(Is[c]),
                                 "CPS": float(cps.get(donor, np.nan))})
        print(f"[{st}] sections used: "
              f"{len(set(r['section'] for r in sec_rows if r['subtype']==st))}")

    sec = pd.DataFrame(sec_rows)
    sec.to_csv(OUT / "sq3_persection_morans.csv", index=False)

    # ---- per-PC correlations (top-5 spatial PCs per subtype) ----
    corr_rows = []
    for st in ac.SUBTYPES:
        pdd = perdonor[perdonor.subtype == st]
        strength = pdd.groupby("pc")["pooled_I"].first().sort_values(ascending=False)
        topk = list(strength.head(TOP_K).index)
        for pc in topk:
            # donor-level (validated)
            dsub = pdd[pdd.pc == pc].dropna(subset=["morans_I", "CPS"])
            if len(dsub) >= 4 and dsub.morans_I.nunique() > 1 and dsub.CPS.nunique() > 1:
                dr, dp = spearmanr(dsub.morans_I, dsub.CPS)
            else:
                dr, dp = np.nan, np.nan
            # section-level
            ssub = sec[(sec.subtype == st) & (sec.pc == pc)].dropna(subset=["morans_I", "CPS"])
            if len(ssub) >= 4 and ssub.morans_I.nunique() > 1 and ssub.CPS.nunique() > 1:
                sr, _ = spearmanr(ssub.morans_I, ssub.CPS)
                dmap = ssub.groupby("donor")["CPS"].first().to_dict()
                sp_p = donor_perm_p(ssub.morans_I.values, ssub.donor.values,
                                    dmap, sr, rng)
            else:
                sr, sp_p = np.nan, np.nan
            corr_rows.append({"subtype": st, "pc": pc,
                              "pooled_I": float(strength.get(pc, np.nan)),
                              "donor_r": dr, "donor_p": dp,
                              "sec_r": sr, "sec_perm_p": sp_p,
                              "n_sections": len(ssub),
                              "n_donors": ssub.donor.nunique() if len(ssub) else 0})
    cdf = pd.DataFrame(corr_rows)
    cdf["sec_fdr_q"] = np.nan
    for st in ac.SUBTYPES:
        m = cdf.subtype == st
        cdf.loc[m, "sec_fdr_q"] = bh_fdr(cdf.loc[m, "sec_perm_p"].values)
    cdf.to_csv(OUT / "sq3_persection_corr.csv", index=False)

    print("\n===== PER-SECTION validation (top-5 spatial PCs) =====")
    print(cdf.round(4).to_string(index=False))

    # ---- overlay scatter plots ----
    for st in ac.SUBTYPES:
        pdd = perdonor[perdonor.subtype == st]
        strength = pdd.groupby("pc")["pooled_I"].first().sort_values(ascending=False)
        topk = list(strength.head(TOP_K).index)
        donors_all = sorted(pdd["donor"].unique())
        cmap = plt.get_cmap("tab10" if len(donors_all) <= 10 else "tab20")
        dcol = {d: cmap(i % cmap.N) for i, d in enumerate(donors_all)}

        fig, axs = plt.subplots(1, len(topk), figsize=(3.6 * len(topk), 3.8), squeeze=False)
        axs = axs[0]
        for ax, pc in zip(axs, topk):
            dsub = pdd[pdd.pc == pc].dropna(subset=["morans_I", "CPS"])
            ssub = sec[(sec.subtype == st) & (sec.pc == pc)].dropna(subset=["morans_I", "CPS"])
            # section points (faded), coloured by donor
            for _, r in ssub.iterrows():
                ax.scatter(r["CPS"], r["morans_I"], s=28, alpha=0.35,
                           color=dcol.get(r["donor"], "0.5"), zorder=2)
            # donor points (solid), coloured by donor
            for _, r in dsub.iterrows():
                ax.scatter(r["CPS"], r["morans_I"], s=95, alpha=1.0,
                           color=dcol.get(r["donor"], "0.5"),
                           edgecolor="black", linewidth=0.6, zorder=4)
            # trendlines
            cr = cdf[(cdf.subtype == st) & (cdf.pc == pc)]
            xs = np.linspace(min(dsub.CPS.min(), ssub.CPS.min()),
                             max(dsub.CPS.max(), ssub.CPS.max()), 50)
            if len(dsub) >= 2:
                a, b = np.polyfit(dsub.CPS, dsub.morans_I, 1)
                ax.plot(xs, a * xs + b, "-", color="black", lw=1.6, zorder=3,
                        label="donor trend")
            if len(ssub) >= 2:
                a2, b2 = np.polyfit(ssub.CPS, ssub.morans_I, 1)
                ax.plot(xs, a2 * xs + b2, "--", color="0.4", lw=1.4, zorder=3,
                        label="section trend")
            ax.axhline(0, color="0.7", lw=0.6)
            t = (f"PC{pc} (pooled I={cr.iloc[0]['pooled_I']:.2f})\n"
                 f"donor r={cr.iloc[0]['donor_r']:.2f} p={cr.iloc[0]['donor_p']:.3f}\n"
                 f"sec r={cr.iloc[0]['sec_r']:.2f} permp={cr.iloc[0]['sec_perm_p']:.3f}"
                 ) if len(cr) else f"PC{pc}"
            ax.set_title(t, fontsize=8)
            ax.set_xlabel("CPS")
        axs[0].set_ylabel("Moran's I")
        axs[0].legend(fontsize=6, loc="best")
        fig.suptitle(f"SQ3 donor (solid) vs section (faded) — {st}  "
                     f"(1 colour/donor)", fontsize=11)
        fig.tight_layout()
        fig.savefig(OUT / f"sq3_overlay_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)
    print("sq3_persection done ->", OUT)


if __name__ == "__main__":
    main()
