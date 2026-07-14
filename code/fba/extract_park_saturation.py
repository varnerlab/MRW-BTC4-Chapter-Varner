"""
Extract the urea-cycle enzyme [substrate]/Km rows from the Park et al. 2016
supplement (Nat. Chem. Biol., PMC4912430), Supplementary Data Set
"Comparison of absolute concentrations to Km"
(data/41589_2016_BFnchembio2077_MOESM585_ESM.xlsx), into a small committed CSV.

Run once (needs openpyxl; the raw supplement is gitignored):
    python3 -m venv .venv && .venv/bin/pip install openpyxl pandas
    .venv/bin/python code/fba/extract_park_saturation.py
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(ROOT, "data", "41589_2016_BFnchembio2077_MOESM585_ESM.xlsx")
OUT = os.path.join(ROOT, "code", "data", "park_saturation.csv")

# EC number -> reaction label in the urea-cycle model
EC2RX = {
    "6.3.4.5": "v1",     # argininosuccinate synthetase
    "4.3.2.1": "v2",     # argininosuccinate lyase (only products present)
    "3.5.3.1": "v3",     # arginase
    "2.1.3.3": "v4",     # ornithine transcarbamylase
    "1.14.13.39": "v5",  # nitric oxide synthase
}

df = pd.read_excel(SRC, header=1)
df.columns = [str(c).strip() for c in df.columns]
df["EC Number"] = df["EC Number"].astype(str).str.strip()
sub = df[df["EC Number"].isin(EC2RX)].copy()
sub["reaction"] = sub["EC Number"].map(EC2RX)
out = sub[["reaction", "EC Number", "Metabolite", "Concentration (M)", "Km (M)*", "Organism"]]
out.columns = ["reaction", "ec", "substrate", "conc_M", "km_M", "organism"]
out = out.sort_values(["reaction", "substrate"]).reset_index(drop=True)

with open(OUT, "w") as f:
    f.write("# Urea-cycle enzyme substrate concentration / Km pairs\n")
    f.write("# Source: Park et al. 2016, Nat. Chem. Biol. 12:482-489 (PMC4912430),\n")
    f.write("#   Supplementary Data Set 'Comparison of absolute concentrations to Km',\n")
    f.write("#   file 41589_2016_BFnchembio2077_MOESM585_ESM.xlsx (concentration + BRENDA Km, in M).\n")
    f.write("# NOTE: EC 4.3.2.1 (v2) lists only the reverse-direction PRODUCTS (arginine, fumarate);\n")
    f.write("#   the forward substrate argininosuccinate is absent from the supplement.\n")
    out.to_csv(f, index=False)

print(f"wrote {OUT} with {len(out)} rows")
