"""
SQ2 / variogram -- the headline figure.

For each subtype, pool all WITHIN-SECTION cell pairs, bin them by spatial
distance, and plot mean GRN dissimilarity (1 - cosine) vs distance. A rising
curve = GRNs diverge with spatial distance. The expression-space curve
(1 - cosine on the 500 HVGs) is overlaid as the circularity baseline: the GRN
curve must rise *beyond* the expression curve to be more than a mapping artefact.

Outputs -> analysis/SQ2/variogram/
  variogram_<slug>.csv          per-bin: dist_center, mean_dissim, sem, n_pairs (GRN + expr)
  variogram_<slug>.png          per-subtype GRN vs expression curve
  variogram_all.png             all subtypes, GRN curves on one axis
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "variogram")
N_BINS = 20


def binned_curve(sdist, dissim, edges):
    """Mean dissimilarity per distance bin."""
    which = np.digitize(sdist, edges) - 1
    centers, means, sems, counts = [], [], [], []
    for b in range(len(edges) - 1):
        m = which == b
        if m.sum() < 10:
            continue
        d = dissim[m]
        centers.append(0.5 * (edges[b] + edges[b + 1]))
        means.append(d.mean())
        sems.append(d.std(ddof=1) / np.sqrt(len(d)))
        counts.append(int(m.sum()))
    return np.array(centers), np.array(means), np.array(sems), np.array(counts)


def main():
    rng = np.random.default_rng(ac.SEED)
    X, cell_ids, meta, expr = ac.load_cache(with_expr=True)
    summary_curves = {}

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) == 0:
            print(f"[{st}] no cells, skipping"); continue

        grn_blocks = ac.section_blocks(X, meta, rows, rng)
        expr_blocks = ac.section_blocks(expr, meta, rows, rng)
        if not grn_blocks:
            print(f"[{st}] no sections pass min-cell cutoff, skipping"); continue

        sd_g, dd_g = ac.pooled_pairs(grn_blocks)
        sd_e, dd_e = ac.pooled_pairs(expr_blocks)

        # quantile bin edges from the GRN spatial distances (comparable counts)
        edges = np.quantile(sd_g, np.linspace(0, 1, N_BINS + 1))
        edges = np.unique(edges)
        cg, mg, eg, ng = binned_curve(sd_g, dd_g, edges)
        ce, me, ee, ne = binned_curve(sd_e, dd_e, edges)

        df = pd.DataFrame({"dist_center": cg, "grn_mean_dissim": mg,
                           "grn_sem": eg, "grn_n_pairs": ng})
        de = pd.DataFrame({"dist_center": ce, "expr_mean_dissim": me,
                           "expr_sem": ee, "expr_n_pairs": ne})
        out = pd.merge(df, de, on="dist_center", how="outer").sort_values("dist_center")
        out.to_csv(OUT / f"variogram_{ac.SLUG[st]}.csv", index=False)
        summary_curves[st] = (cg, mg)

        # per-subtype figure: GRN vs expression baseline
        fig, axL = plt.subplots(figsize=(6, 4.5))
        axL.errorbar(cg, mg, yerr=eg, fmt="-o", ms=3, lw=1.5,
                     color="C0", label="GRN (1 - cosine)")
        axL.errorbar(ce, me, yerr=ee, fmt="--s", ms=3, lw=1.2,
                     color="C3", label="expression baseline")
        axL.set_xlabel("within-section spatial distance (µm)")
        axL.set_ylabel("mean dissimilarity (1 - cosine)")
        axL.set_title(f"SQ2 variogram — {st}\n"
                      f"{len(rows):,} cells, {len(grn_blocks)} sections, "
                      f"{len(sd_g):,} pairs")
        axL.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(OUT / f"variogram_{ac.SLUG[st]}.png", dpi=150)
        plt.close(fig)
        print(f"[{st}] {len(rows):,} cells | {len(grn_blocks)} sections | "
              f"{len(sd_g):,} pairs | GRN dissim {mg.min():.3f}->{mg.max():.3f}")

    # combined GRN figure
    if summary_curves:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for i, (st, (c, m)) in enumerate(summary_curves.items()):
            ax.plot(c, m, "-o", ms=3, lw=1.5, color=f"C{i}", label=st)
        ax.set_xlabel("within-section spatial distance (µm)")
        ax.set_ylabel("mean GRN dissimilarity (1 - cosine)")
        ax.set_title("SQ2 variogram — all subtypes")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(OUT / "variogram_all.png", dpi=150)
        plt.close(fig)
    print("variogram done ->", OUT)


if __name__ == "__main__":
    main()
