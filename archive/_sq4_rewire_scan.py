"""Phase 1: rewiring scan over the top-20 robust L2/3 regulators.

For each gene, between SQ2 spatial clusters 0/1, measure whether the regulon's SHAPE
differs (rewiring) vs just its loudness (on/off). Shape = L1-normalised mean regulon
over ACTIVE cells in each cluster; cosine of the two = how similar the targeting is.
Low cosine + enough active cells in BOTH clusters = a genuine rewiring candidate.
"""
import numpy as np
import pandas as pd
import analysis_common as ac

ST = "L2/3 IT"
TOPN = 20
MIN_ACT = 100   # active cells needed per cluster for a trustworthy shape
OUT = ac.ANALYSIS / "SQ4" / "sq4_rewire_scan.csv"
LAB_CSV = ac.ANALYSIS / "SQ2" / "grn_clustering" / "clustering_labels_l23_it.csv"

X, cell_ids, meta, _ = ac.load_cache()
idx = np.load(ac.CACHE / "spatial_grn_index.npz", allow_pickle=True)
edge_src = idx["edge_src"].astype(str)

master = pd.read_csv(ac.ANALYSIS / "SQ4" / "sq4_regulator_morans.csv")
m = master[master.subtype == ST]
rob = m[(m.q_fdr < 0.05) & (m.q_eqsec < 0.05)].sort_values("I_R", ascending=False)
top_genes = rob.head(TOPN)["regulator"].tolist()
print(f"top {len(top_genes)} robust L2/3 genes: {top_genes}")

rows = ac.subtype_rows(meta, ST)
cids = cell_ids[rows]
lab_df = pd.read_csv(LAB_CSV)
lab_map = dict(zip(lab_df.cell_id.astype(str), lab_df.grn_cluster.astype(int)))
lab = np.array([lab_map.get(c, -1) for c in cids])
valid = lab >= 0

def l1norm(v):
    s = np.abs(v).sum()
    return v / s if s > 0 else v

recs = []
for R in top_genes:
    cols = np.flatnonzero(edge_src == R)
    M = X[rows][:, cols].toarray().astype(np.float64)
    irow = m[m.regulator == R].iloc[0]
    rec = {"regulator": R, "I_R": float(irow.I_R),
           "q_fdr": float(irow.q_fdr), "q_eqsec": float(irow.q_eqsec)}
    cl = {}
    for c in (0, 1):
        sub = M[(lab == c) & valid]
        act = np.abs(sub).sum(1) > 0
        cl[c] = {"nact": int(act.sum()),
                 "mean_active": sub[act].mean(0) if act.sum() > 0 else None,
                 "mag": float(np.abs(sub.mean(0)).sum())}
    rec["n_active_c0"], rec["n_active_c1"] = cl[0]["nact"], cl[1]["nact"]
    rec["mag_c0"], rec["mag_c1"] = cl[0]["mag"], cl[1]["mag"]
    hi, lo = max(cl[0]["mag"], cl[1]["mag"]), min(cl[0]["mag"], cl[1]["mag"])
    rec["mag_ratio"] = hi / (lo + 1e-12)
    if cl[0]["nact"] >= MIN_ACT and cl[1]["nact"] >= MIN_ACT:
        s0, s1 = l1norm(cl[0]["mean_active"]), l1norm(cl[1]["mean_active"])
        cos = float(s0 @ s1 / (np.linalg.norm(s0) * np.linalg.norm(s1) + 1e-12))
        rec["shape_cosine"] = cos
        rec["rewiring_score"] = 1.0 - cos
        rec["both_active"] = True
    else:
        rec["shape_cosine"] = np.nan
        rec["rewiring_score"] = np.nan
        rec["both_active"] = False
    recs.append(rec)
    print(f"  {R:<13} I_R={rec['I_R']:.3f} nact=({rec['n_active_c0']},"
          f"{rec['n_active_c1']}) mag_ratio={rec['mag_ratio']:.2f} "
          f"cos={rec['shape_cosine'] if rec['both_active'] else 'NA'}")

df = pd.DataFrame(recs).sort_values("rewiring_score", ascending=False, na_position="last")
df.to_csv(OUT, index=False)
print("\n=== ranked by rewiring (1 - shape cosine; high = more rewiring) ===")
print(df[["regulator", "I_R", "n_active_c0", "n_active_c1", "mag_ratio",
          "shape_cosine", "rewiring_score", "both_active"]].to_string(index=False))
print("saved", OUT)
