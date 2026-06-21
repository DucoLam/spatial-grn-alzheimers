import pandas as pd
import analysis_common as ac

m = pd.read_parquet(ac.CACHE / "spatial_meta.parquet")
print("COLS:", list(m.columns))
print("NROWS:", len(m))
print("N unique sections:", m["section"].astype(str).nunique())

# find a donor-like column
cand = [c for c in m.columns if "donor" in c.lower() or "specimen" in c.lower()
        or "subject" in c.lower() or "individual" in c.lower()]
print("donor-like cols:", cand)

dcol = cand[0] if cand else None
if dcol is not None:
    print(f"\nN unique donors ({dcol}):", m[dcol].nunique())
    g = m.groupby(dcol)["section"].nunique().sort_values(ascending=False)
    print(f"\nsections per donor:\n{g.to_string()}")
    print("\nTOTAL sections:", int(g.sum()), "| donors:", g.shape[0])
else:
    # fall back: maybe section name encodes donor; show all sections
    secs = sorted(m["section"].astype(str).unique())
    print("\nall sections:")
    for s in secs:
        print("  ", s)
