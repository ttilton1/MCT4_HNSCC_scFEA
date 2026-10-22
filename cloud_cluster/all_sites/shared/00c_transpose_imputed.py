#!/usr/bin/env python
# ==============================================================================
# 00c_transpose_imputed.py  --  make the shared imputed matrices R-readable.
#
# 00b writes imputed_expression.csv as genes x cells (661 rows x up-to-244K
# cols). R's data.table::fread mis-parses that width and returns garbage numbers
# (Python/pandas reads it fine, so scFEA predict is unaffected). This script
# reads each wide matrix in pandas (correct) and writes a TALL copy,
# imputed_tall.csv = cells x genes (up-to-244K rows x 662 cols: barcode + 661
# genes), which fread reads reliably. Both bulk_before and bulk_after 01 scripts
# read the tall file.
#
# Keeps the wide file (scFEA predict for bulk_after reads it directly in Python).
# Values are unchanged -- this is a pure reshape. Prints value range per type as
# a sanity check (expect ~0 to ~5 for MAGIC-imputed log-normalized data).
# ==============================================================================
import os
import sys
import numpy as np
import pandas as pd

SHARED_BASE = "/project2/dechenli_1477/ttilton/spring26/scfea/data/input/all_sites/shared"
TYPES = ["T","Myeloid","Fibroblast","Epithelial","Endothelial","B","NK","Mast","Plasma"]

def main():
    for ct in TYPES:
        wide = os.path.join(SHARED_BASE, ct, "imputed_expression.csv")
        tall = os.path.join(SHARED_BASE, ct, "imputed_tall.csv")
        if not os.path.exists(wide):
            print(f"{ct}: MISSING {wide}, skipping", flush=True); continue

        # genes x cells (index=genes, columns=barcodes)
        df = pd.read_csv(wide, index_col=0)
        v = df.values
        lo, hi = float(np.nanmin(v)), float(np.nanmax(v))
        bad = bool(np.isnan(v).any() or np.isinf(v).any() or lo < -1e-6 or hi > 100)
        print(f"{ct}: {df.shape[0]} genes x {df.shape[1]} cells | "
              f"range [{lo:.4f}, {hi:.4f}] | {'SUSPECT' if bad else 'ok'}", flush=True)
        if bad:
            print(f"  WARNING: {ct} values look off in the SOURCE wide file -- check 00b", flush=True)

        # transpose -> cells x genes, write tall (barcode index, gene columns)
        df.T.to_csv(tall)  # index label defaults to the barcodes' axis name (blank) -> fine
        gb = os.path.getsize(tall) / (1024**3)
        print(f"  wrote {tall} ({df.shape[1]} cells x {df.shape[0]} genes, {gb:.2f} GB)", flush=True)
        del df

    print("done.", flush=True)

if __name__ == "__main__":
    main()

