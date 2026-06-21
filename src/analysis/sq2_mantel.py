"""
SQ2 / mantel -- the headline statistic.

Stratified (within-section) Mantel test per subtype: correlate the pooled
within-section spatial-distance vector with the GRN-dissimilarity vector.
Significance via restricted permutation (cell labels shuffled within each
section, preserving the per-section coordinate frames). The same test is run
on the expression baseline so the GRN signal can be compared against what raw
expression already explains.

Positive r with small p  =>  cells farther apart spatially have more different
GRNs (spatial structure in GRNs).

Outputs -> analysis/SQ2/mantel/
  mantel_results.csv     one row per subtype x feature-space (GRN, expression)
  mantel_persection.csv  per-section Pearson r (descriptive)
"""
import numpy as np
import pandas as pd

import analysis_common as ac

OUT = ac.ensure_dir(ac.ANALYSIS / "SQ2" / "mantel")


def per_section_r(blocks):
    rows = []
    for b in blocks:
        iu = np.triu_indices(b["n"], k=1)
        r = ac._pearson(b["S"][iu], b["D"][iu])
        rows.append({"section": b["section"], "n_cells": b["n"],
                     "n_pairs": len(iu[0]), "r_pearson": r})
    return rows


def main():
    X, cell_ids, meta, expr = ac.load_cache(with_expr=True)
    results, persection = [], []

    for st in ac.SUBTYPES:
        rows = ac.subtype_rows(meta, st)
        if len(rows) == 0:
            print(f"[{st}] no cells, skipping"); continue

        for space, mat in (("GRN", X), ("expression", expr)):
            rng = np.random.default_rng(ac.SEED)   # same subsample per space
            blocks = ac.section_blocks(mat, meta, rows, rng)
            if not blocks:
                print(f"[{st}/{space}] no sections pass cutoff, skipping"); continue
            res = ac.stratified_mantel(blocks, rng)
            res.update({"subtype": st, "space": space, "n_cells": int(len(rows))})
            results.append(res)
            print(f"[{st}/{space}] r={res['r_pearson']:.4f} "
                  f"rho={res['r_spearman']:.4f} p={res['p_perm']:.4f} "
                  f"({res['n_pairs']:,} pairs, {res['n_sections']} sections)")
            if space == "GRN":
                for r in per_section_r(blocks):
                    r["subtype"] = st
                    persection.append(r)

    cols = ["subtype", "space", "n_cells", "n_sections", "n_pairs",
            "r_pearson", "r_spearman", "p_perm"]
    pd.DataFrame(results)[cols].to_csv(OUT / "mantel_results.csv", index=False)
    pd.DataFrame(persection).to_csv(OUT / "mantel_persection.csv", index=False)
    print("mantel done ->", OUT)


if __name__ == "__main__":
    main()
