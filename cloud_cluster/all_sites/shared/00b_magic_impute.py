#!/usr/bin/env python
# ==============================================================================
# 00b_magic_impute.py  --  run MAGIC ONCE on the full normalized cell matrix
#                          for one major type. Shared by Run A and Run B.
#
# Reads (from 00a, per major type):
#   normalized.mtx   genes x cells (log-normalized)
#   genes.txt        gene names (mtx row order)
#   barcodes.txt     cell barcodes (mtx column order)
#
# Runs magic.MAGIC(solver=...).fit_transform on cells x genes. The graph is
# built on ALL genes (matches scFEA's internal MAGIC), but only the ~150 scFEA
# metabolic OVERLAP genes are RETURNED -- that is all scFEA uses, and it keeps
# the output tiny even for T (244K cells). Returning a subset does not change
# the imputed values: the diffusion operator is built from the full input
# regardless; `genes` only selects which columns the operator is applied to.
#
# Writes:
#   imputed_expression.csv   overlap genes x cells (scFEA input for BOTH runs)
#
# scFEA is then run downstream with --sc_imputation False (imputation done here,
# once). To match scFEA's internal call exactly, every MAGIC parameter is left
# at magic.MAGIC() defaults (knn=5, decay=1, t='auto', n_pca=100) EXCEPT solver,
# which is exposed: 'exact' matches scFEA but is memory-heavy at ~244K cells;
# 'approximate' is MAGIC's recommended large-data setting. Pick one solver for
# ALL major types so imputation is uniform across types.
#
# Usage:
#   python 00b_magic_impute.py \
#       --combo_dir  /.../data/input/all_sites_imputed/Myeloid \
#       --module_file /.../scFEA_src/data/module_gene_m168.csv \
#       --solver exact
# ==============================================================================
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import magic


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combo_dir", required=True,
                    help="dir with normalized.mtx, genes.txt, barcodes.txt")
    ap.add_argument("--module_file", required=True,
                    help="scFEA module_gene_m168.csv (defines overlap genes)")
    ap.add_argument("--solver", default="exact", choices=["exact", "approximate"],
                    help="MAGIC solver. exact matches scFEA; approximate scales to large N.")
    ap.add_argument("--out", default=None,
                    help="output CSV (default: <combo_dir>/imputed_expression.csv)")
    args = ap.parse_args()

    combo = os.path.basename(os.path.normpath(args.combo_dir))
    out_path = args.out or os.path.join(args.combo_dir, "imputed_expression.csv")

    print(f"=== MAGIC once | combo={combo} | solver={args.solver} ===", flush=True)

    # ---- load sparse normalized matrix (genes x cells) ----
    print("Reading normalized.mtx ...", flush=True)
    mat = scipy.io.mmread(os.path.join(args.combo_dir, "normalized.mtx")).tocsr()  # genes x cells
    with open(os.path.join(args.combo_dir, "genes.txt")) as f:
        genes = [ln.strip() for ln in f]
    with open(os.path.join(args.combo_dir, "barcodes.txt")) as f:
        barcodes = [ln.strip() for ln in f]
    assert mat.shape == (len(genes), len(barcodes)), \
        f"mtx shape {mat.shape} != ({len(genes)} genes, {len(barcodes)} cells)"
    print(f"  {mat.shape[0]} genes x {mat.shape[1]} cells", flush=True)

    # ---- MAGIC convention is cells x genes ----
    X = mat.T.tocsr()  # cells x genes

    # ---- overlap genes: which columns MAGIC should return ----
    module_df = pd.read_csv(args.module_file, index_col=0)
    scfea_genes = set(
        str(g) for g in module_df.values.ravel()
        if pd.notna(g) and str(g) != "" and str(g) != "nan"
    )
    gene_pos = {g: i for i, g in enumerate(genes)}
    overlap = [g for g in genes if g in scfea_genes]         # preserves mtx order
    overlap_idx = [gene_pos[g] for g in overlap]
    print(f"  overlap metabolic genes to return: {len(overlap)}", flush=True)
    if len(overlap) == 0:
        sys.exit("No overlap between module genes and matrix genes -- check gene naming.")

    # ---- MAGIC: defaults except solver (matches scFEA's magic.MAGIC()) ----
    t0 = time.time()
    op = magic.MAGIC(solver=args.solver)
    print(f"Running MAGIC.fit_transform (solver={args.solver}) ...", flush=True)
    imputed = op.fit_transform(X, genes=overlap_idx)          # cells x len(overlap), dense
    dt = time.time() - t0
    print(f"  MAGIC done in {dt/60:.1f} min, imputed shape {np.asarray(imputed).shape}", flush=True)

    # ---- assemble genes x cells and write ----
    imp = np.asarray(imputed)
    df = pd.DataFrame(imp, index=barcodes, columns=overlap).T  # genes x cells
    df.to_csv(out_path)
    gb = os.path.getsize(out_path) / (1024 ** 3)
    print(f"  Wrote {out_path} ({df.shape[0]} genes x {df.shape[1]} cells, {gb:.3f} GB)", flush=True)
    print("=== done ===", flush=True)


if __name__ == "__main__":
    main()

