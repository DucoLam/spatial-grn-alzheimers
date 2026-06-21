"""
SQ3 (revised) -- per-donor spatial GRN organisation vs AD severity.

For each donor x subtype: compute the degree of spatial GRN organisation
as the MEDIAN Moran's I across the top N_PCS GRN principal components
(using the within-section kNN spatial-weights graph, same setup as SQ2).
Then correlate that per-donor score with each severity axis (Braak, Thal,
CERAD) using Spearman across donors. BH-FDR across the focused test set.

This directly answers SQ3: "do donors with higher AD severity show stronger
(or weaker) spatial GRN clustering in their MTG cells?"

No permutation test per-donor (would be 9 x 4 x 999 = slow); instead we
use the observed Moran's I per PC as a summary statistic. The N=donors
Spearman is the outer (cross-donor) test.

Outputs -> analysis/SQ3/
  sq3_spatial_severity_scores.csv   donor x subtype: median_I, mean_I, n_sig_pcs, severity cols
  sq3_spatial_severity_corr.csv     subtype x severity: spearman_r, p, fdr_q, N
  sq3_spatial_severity_<slug>.png   scatter: median_I vs each severity axis
"""
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import scipy.sparse as sp
from scipy.stats import spearmanr

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ3")
N_PCS   = ac.N_PCS          # 20
KNN_K   = ac.KNN_K          # 6
MIN_CELLS_PER_SECTION = ac.MIN_CELLS_PER_SECTION   # 30
MIN_DONOR_CELLS = 50        # minimum cells per donor x subtype to include

SEVERITY_AXES = ["CPS"]

CERAD_MAP = {"Absent": 0, "Sparse": 1, "Moderate": 2, "Frequent": 3}
def braak_n(v): return {"0":0,"I":1,"II":2,"III":3,"IV":4,"V":5,"VI":6}.get(str(v).replace("Braak","").strip(), np.nan)
def thal_n(v):
    try: return float(str(v).replace("Thal","").strip())
    except Exception: return np.nan


def bh_fdr(pvals):
    p = np.asarray(pvals, float); ok = ~np.isnan(p)
    q = np.full(len(p), np.nan)
    pv = p[ok]; n = len(pv)
    if n == 0: return q
    order = np.argsort(pv)
    ranked = np.empty(n); ranked[order] = np.arange(1, n + 1)
    qv = pv * n / ranked
    qs = qv[order]
    qs = np.minimum.accumulate(qs[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(qs, 0, 1)
    q[ok] = out
    return q


def donor_severity():
    rna = ad.read_h5ad(
        "/tudelft.net/staff-umbrella/ScReNI/dflam/data/seaad_paired_rna_pool4.h5ad",
        backed="r")
    obs = rna.obs.copy()
    obs["Donor ID"] = obs["Donor ID"].astype(str)
    g = obs.groupby("Donor ID", observed=True).first()
    s = pd.DataFrame(index=g.index)
    s["Braak"] = g["Braak"].map(braak_n)
    s["Thal"]  = g["Thal"].map(thal_n)
    s["CERAD"] = g["CERAD score"].astype(str).map(CERAD_MAP)
    s["CPS"]   = pd.to_numeric(g["Continuous Pseudo-progression Score"], errors="coerce")
    s["MMSE"]  = pd.to_numeric(g["Last MMSE Score"], errors="coerce")
    s = s.apply(lambda col: pd.to_numeric(col, errors="coerce"))
    return s


def within_section_knn_W_donor(coords, sections, k=KNN_K):
    """Row-standardised kNN W for a subset of cells (one donor x subtype)."""
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
    data = np.ones(len(rptr), dtype=np.float64)
    W = sp.csr_matrix((data, (rptr, cidx)), shape=(n, n))
    rs = np.asarray(W.sum(1)).ravel()
    rs[rs == 0] = 1.0
    W = sp.diags(1.0 / rs) @ W
    return W


def donor_morans_I(G, W):
    """
    Moran's I for each PC in G (n_cells x n_pcs) using weights W.
    Returns array of shape (n_pcs,) — NaN if W is None or degenerate.
    """
    n_pcs = G.shape[1]
    Is = np.full(n_pcs, np.nan)
    if W is None:
        return Is
    Wsum = W.sum()
    n = G.shape[0]
    if Wsum == 0:
        return Is
    for c in range(n_pcs):
        z = G[:, c] - G[:, c].mean()
        den = float(z @ z)
        if den == 0:
            continue
        num = float(z @ (W @ z))
        Is[c] = (n / Wsum) * (num / den)
    return Is


def load_sq2_pooled_morans(slug, n_comp):
    """
    Pooled (severity-blind) GRN Moran's I per PC from the SQ2 output, aligned to
    PC index 0..n_comp-1. Returns (pooled_I[n_comp], sig_mask[n_comp]).
    PCs match because both analyses run TruncatedSVD(seed) on the same subtype X.
    """
    path = ac.ANALYSIS / "SQ2" / "morans_i" / f"morans_{slug}.csv"
    pooled_I = np.full(n_comp, np.nan)
    sig_mask = np.zeros(n_comp, dtype=bool)
    if not path.exists():
        return pooled_I, sig_mask
    d = pd.read_csv(path)
    d = d[d["space"] == "GRN"]
    for _, r in d.iterrows():
        i = int(r["pc"]) - 1          # SQ2 pc is 1-based
        if 0 <= i < n_comp:
            pooled_I[i] = float(r["morans_I"])
            sig_mask[i] = bool(r["p_perm"] < 0.05)
    return pooled_I, sig_mask


def main():
    X, cell_ids, meta, _ = ac.load_cache(with_expr=False)
    meta = meta.copy()
    meta["Donor ID"] = meta["Donor ID"].astype(str)
    sev = donor_severity()

    score_rows = []
    perpc_rows = []          # long format: one row per (subtype, donor, PC)
    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        meta_st = meta.iloc[rows].reset_index(drop=True)
        X_st = X[rows]

        # GRN PCs across ALL cells of this subtype (same as SQ2 — global PCs)
        n_comp = min(N_PCS, len(rows) - 1)
        svd = TruncatedSVD(n_components=n_comp, random_state=ac.SEED)
        G_all = svd.fit_transform(X_st)   # shape: (n_cells_st, n_comp)

        # pooled, severity-blind spatial strength per PC (from SQ2)
        pooled_I, sig_mask = load_sq2_pooled_morans(ac.SLUG[st], n_comp)
        # rank PCs by pooled spatial strength (descending); ignore NaNs
        order = np.argsort(np.where(np.isnan(pooled_I), -np.inf, pooled_I))[::-1]
        topI5_idx  = order[:5]
        topI10_idx = order[:10]
        w_pooled = np.clip(np.nan_to_num(pooled_I, nan=0.0), 0, None)  # weights

        donors = meta_st["Donor ID"].values
        for donor in pd.unique(donors):
            mask = donors == donor
            if mask.sum() < MIN_DONOR_CELLS:
                continue
            G_d = G_all[mask]
            coords_d = meta_st[["x", "y"]].values[mask]
            sects_d  = meta_st["section"].values[mask]

            # drop sections with too few cells
            valid_sects = [s for s in pd.unique(sects_d)
                           if (sects_d == s).sum() >= MIN_CELLS_PER_SECTION]
            if not valid_sects:
                continue
            keep = np.isin(sects_d, valid_sects)
            G_d      = G_d[keep]
            coords_d = coords_d[keep]
            sects_d  = sects_d[keep]

            if G_d.shape[0] < MIN_DONOR_CELLS:
                continue

            W_d = within_section_knn_W_donor(coords_d, sects_d)
            Is  = donor_morans_I(G_d, W_d)   # (n_comp,) Moran's I per PC, in PC (variance) order
            var = svd.explained_variance_ratio_[:len(Is)]

            # summaries. Is is ordered by PC index = variance order, so the
            # "top-k" leading PCs are simply the first k entries (NOT selected
            # on Moran's I, which would be circular).
            def topk_mean(k):
                vals = Is[:k]
                return float(np.nanmean(vals)) if np.any(~np.isnan(vals)) else np.nan
            ok = ~np.isnan(Is)
            vw_mean = (float(np.sum(Is[ok] * var[ok]) / np.sum(var[ok]))
                       if np.any(ok) and np.sum(var[ok]) > 0 else np.nan)

            # severity-blind spatial-strength-weighted mean (weights from SQ2)
            wok = ok & (w_pooled > 0)
            wI_mean = (float(np.sum(Is[wok] * w_pooled[wok]) / np.sum(w_pooled[wok]))
                       if np.any(wok) else np.nan)
            def idx_mean(idx):
                vals = Is[idx]
                return float(np.nanmean(vals)) if np.any(~np.isnan(vals)) else np.nan
            mean_sigpcs = (float(np.nanmean(Is[sig_mask]))
                           if sig_mask.any() and np.any(~np.isnan(Is[sig_mask]))
                           else np.nan)

            row = {
                "subtype":    st,
                "donor":      donor,
                "n_cells":    int(G_d.shape[0]),
                "n_sections": len(valid_sects),
                "n_pcs":      int(np.sum(ok)),
                # severity-blind spatial-strength selection (the intended metrics)
                "wI_mean":    wI_mean,
                "mean_topI5": idx_mean(topI5_idx),
                "mean_topI10":idx_mean(topI10_idx),
                "mean_sigpcs":mean_sigpcs,
                # variance-based + plain summaries (for comparison)
                "vw_mean":    vw_mean,
                "mean_top3":  topk_mean(3),
                "mean_top5":  topk_mean(5),
                "mean_top10": topk_mean(10),
                "median_I":   float(np.nanmedian(Is)),
                "mean_I":     float(np.nanmean(Is)),
                "max_I":      float(np.nanmax(Is)),
            }
            for ax in SEVERITY_AXES:
                row[ax] = sev[ax].get(donor, np.nan)
            score_rows.append(row)

            # per-PC records (the intended SQ3 analysis): keep each PC's own
            # per-donor Moran's I, plus its pooled severity-blind spatial strength
            cps_val = row.get("CPS", np.nan)
            for c in range(len(Is)):
                perpc_rows.append({
                    "subtype":  st,
                    "donor":    donor,
                    "pc":       c + 1,
                    "morans_I": float(Is[c]) if not np.isnan(Is[c]) else np.nan,
                    "pooled_I": float(pooled_I[c]) if not np.isnan(pooled_I[c]) else np.nan,
                    "sq2_sig":  bool(sig_mask[c]),
                    "CPS":      cps_val,
                })
            print(f"  [{st}] {donor}: n={G_d.shape[0]}, "
                  f"wI_mean={row['wI_mean']:.3f}, CPS={cps_val}")

    scores = pd.DataFrame(score_rows)
    scores.to_csv(OUT / "sq3_spatial_severity_scores.csv", index=False)

    # --- Spearman correlations: each summary metric vs each severity axis ---
    SUMMARY_METRICS = ["wI_mean", "mean_topI5", "mean_topI10", "mean_sigpcs",
                       "vw_mean", "mean_top3", "mean_top5", "mean_top10",
                       "median_I", "mean_I", "max_I"]
    corr_rows = []
    for metric in SUMMARY_METRICS:
        for st in ac.SUBTYPES:
            df = scores[scores.subtype == st].copy()
            for ax in SEVERITY_AXES:
                sub = df.dropna(subset=[ax, metric])
                if len(sub) >= 4 and sub[ax].nunique() > 1 and sub[metric].nunique() > 1:
                    r, p = spearmanr(sub[metric], sub[ax])
                else:
                    r, p = np.nan, np.nan
                corr_rows.append({"metric": metric, "subtype": st, "severity": ax,
                                   "spearman_r": r, "p": p, "N": len(sub)})

    cdf = pd.DataFrame(corr_rows)
    # BH-FDR computed WITHIN each metric (each metric is its own analysis;
    # don't pool across metrics or the multiple-testing burden is artificial)
    cdf["fdr_q"] = np.nan
    for metric in SUMMARY_METRICS:
        m = cdf.metric == metric
        cdf.loc[m, "fdr_q"] = bh_fdr(cdf.loc[m, "p"].values)
    cdf.to_csv(OUT / "sq3_spatial_severity_corr.csv", index=False)

    print("\n--- SQ3 spatial severity correlations (by metric) ---")
    for metric in SUMMARY_METRICS:
        sub = cdf[cdf.metric == metric].sort_values("p")
        print(f"\n=== {metric} ===")
        print(sub[["subtype", "severity", "spearman_r", "p", "N", "fdr_q"]]
              .to_string(index=False))

    # ================================================================
    # PER-PC ANALYSIS (the intended SQ3): for each spatially-strong PC,
    # does that PC's per-donor Moran's I track CPS across donors?
    #   - PC selection is by POOLED (SQ2) Moran's I = severity-blind.
    #   - one Spearman per PC (NO averaging across PCs).
    #   - FDR within the top-K spatial PCs per subtype.
    # ================================================================
    TOP_K = 5
    pp = pd.DataFrame(perpc_rows)
    pp.to_csv(OUT / "sq3_perpc_morans.csv", index=False)

    perpc_corr = []
    for st in ac.SUBTYPES:
        d = pp[pp.subtype == st]
        if d.empty:
            continue
        # rank PCs by pooled (severity-blind) spatial strength
        strength = d.groupby("pc")["pooled_I"].first().sort_values(ascending=False)
        rank_of = {pc: i + 1 for i, pc in enumerate(strength.index)}
        topk_pcs = set(strength.head(TOP_K).index)
        for pc in sorted(d["pc"].unique()):
            sub = d[d["pc"] == pc].dropna(subset=["morans_I", "CPS"])
            if len(sub) >= 4 and sub["morans_I"].nunique() > 1 and sub["CPS"].nunique() > 1:
                r, p = spearmanr(sub["morans_I"], sub["CPS"])
            else:
                r, p = np.nan, np.nan
            perpc_corr.append({
                "subtype": st, "pc": pc,
                "pooled_I": float(strength.get(pc, np.nan)),
                "spatial_rank": rank_of.get(pc, np.nan),
                "in_topk": pc in topk_pcs,
                "spearman_r_vs_CPS": r, "p": p, "N": len(sub),
            })
    ppc = pd.DataFrame(perpc_corr)
    # FDR within each subtype across the pre-specified top-K spatial PCs only
    ppc["fdr_q"] = np.nan
    for st in ac.SUBTYPES:
        m = (ppc.subtype == st) & (ppc.in_topk)
        if m.any():
            ppc.loc[m, "fdr_q"] = bh_fdr(ppc.loc[m, "p"].values)
    ppc.to_csv(OUT / "sq3_perpc_corr.csv", index=False)

    print("\n===== PER-PC: each spatial PC's per-donor Moran's I vs CPS =====")
    for st in ac.SUBTYPES:
        sub = ppc[(ppc.subtype == st) & (ppc.in_topk)].sort_values("p")
        if sub.empty:
            continue
        print(f"\n=== {st} (top-{TOP_K} spatial PCs by pooled Moran's I) ===")
        print(sub[["pc", "pooled_I", "spearman_r_vs_CPS", "p", "N", "fdr_q"]]
              .round(4).to_string(index=False))

    # per-PC bar plots: Spearman(Moran's I, CPS) per PC, top-K highlighted
    for st in ac.SUBTYPES:
        sub = ppc[ppc.subtype == st].sort_values("pc")
        if sub.empty:
            continue
        fig, ax = plt.subplots(figsize=(8, 4))
        colors = ["#c0392b" if tk else "0.8" for tk in sub["in_topk"]]
        ax.bar(sub["pc"], sub["spearman_r_vs_CPS"], color=colors)
        for _, rr in sub.iterrows():
            if rr["in_topk"] and not np.isnan(rr["p"]) and rr["p"] < 0.05:
                y = rr["spearman_r_vs_CPS"]
                ax.text(rr["pc"], y + (0.03 if y >= 0 else -0.06), "*",
                        ha="center", fontsize=12)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylim(-1.05, 1.05)
        ax.set_xlabel("GRN principal component")
        ax.set_ylabel("Spearman r  (per-donor Moran's I vs CPS)")
        ax.set_title(f"SQ3 per-PC — {st}  (red = top-{TOP_K} spatial PCs, * = p<0.05)")
        fig.tight_layout()
        fig.savefig(OUT / f"sq3_perpc_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)

    # --- scatter plots: wI_mean (aggregate, kept for reference) vs severity ---
    PLOT_METRIC = "wI_mean"
    for st in ac.SUBTYPES:
        df = scores[scores.subtype == st].copy()
        axes_with_data = [ax for ax in SEVERITY_AXES if df[ax].notna().any()]
        if not axes_with_data:
            continue
        fig, axs = plt.subplots(1, len(axes_with_data),
                                figsize=(4 * len(axes_with_data), 4))
        if len(axes_with_data) == 1:
            axs = [axs]
        for ax_plot, ax_name in zip(axs, axes_with_data):
            sub = df.dropna(subset=[ax_name, PLOT_METRIC])
            ax_plot.scatter(sub[ax_name], sub[PLOT_METRIC], s=60, zorder=3)
            for _, row in sub.iterrows():
                ax_plot.annotate(row["donor"][-6:],
                                 (row[ax_name], row[PLOT_METRIC]),
                                 fontsize=6, xytext=(3, 3),
                                 textcoords="offset points")
            row_c = cdf[(cdf.metric == PLOT_METRIC) & (cdf.subtype == st)
                        & (cdf.severity == ax_name)]
            if len(row_c):
                r = row_c.iloc[0]["spearman_r"]; p = row_c.iloc[0]["p"]
                q = row_c.iloc[0]["fdr_q"]
                ax_plot.set_title(f"{ax_name}\nr={r:.2f} p={p:.3f} q={q:.2f}",
                                  fontsize=9)
            ax_plot.set_xlabel(ax_name)
            ax_plot.set_ylabel("Spatial-weighted mean Moran's I (GRN)")
        fig.suptitle(f"SQ3 spatial organisation vs severity — {st}", fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / f"sq3_spatial_severity_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)

    print("sq3_spatial_severity done ->", OUT)


if __name__ == "__main__":
    main()
