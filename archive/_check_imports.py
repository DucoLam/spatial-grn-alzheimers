import importlib
mods = ["numpy", "scipy", "pandas", "sklearn", "matplotlib", "anndata"]
for m in mods:
    try:
        mod = importlib.import_module(m)
        print("OK", m, getattr(mod, "__version__", "?"))
    except Exception as e:
        print("MISSING", m, repr(e))
