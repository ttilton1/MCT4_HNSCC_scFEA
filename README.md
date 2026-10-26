# MCT4-associated lactate export in HNSCC myeloid cells

Analysis code for:

> Tilton T, Park YM, Zhu H, Zhang C, Lin D-C. *MCT4-Associated Lactate Export Defines
> Myeloid Metabolic Reprogramming during Head and Neck Cancer Progression.*
> bioRxiv 2026. https://doi.org/10.64898/2026.10.13.758108

We infer metabolic flux with [scFEA](https://github.com/changwn/scFEA) across 168 modules
(M1-M168) of the M171 network, in 645,351 cells from 114 patients across five published
HNSCC cohorts. The analysis is built on **metabolic allocation**, defined as a module's
flux as a fraction of its unit's total flux, so we can compare across cell types. Raw flux
will not work for that, because each major cell type has its own fitted model.

---

## Setup

Scripts find the data directory through an environment variable:

```bash
export SCFEA_ROOT=/path/to/your/scfea_data
```

To check a rebuild against ours, run `bash check_figure_scripts.sh`. It verifies that
every input exists, that the matched-null patch ran, and that the figures and tables
regenerate. It does not recompute flux.

Expected layout under `$SCFEA_ROOT`:

```
data/combined_wide/all_sites/<arm>/    assembled flux tables
data/input/all_sites/<arm>/<Type>/     train_expression.csv, per cell type
data/independent/                      GSE139324_RAW/, GSE164690_RAW/
metadata/all_sites/<arm>/<Type>/       metadata_units.csv
output/all_sites/<arm>/<Type>/         raw scFEA flux
results/                               everything written by these scripts
tcga/                                  cBioPortal HNSC PanCancer Atlas download
```

**Requirements.** R 4.2 or later with dplyr, tidyr, readr, data.table, Matrix, ggplot2,
patchwork, ggtext, survival, edgeR, Seurat. Python 3.7 with `magic-impute` and torch for
the cluster steps. A separate environment with `cptac==1.5.14` for the proteomics. See
`envs/` for what we actually ran in.

**Data not included here.** The five source cohorts (He, Kim, Koh, Pittet, Tirosh) as a
harmonized Seurat object. TCGA-HNSC from cBioPortal. CPTAC-HNSCC through the `cptac`
package. GSE139324 (Cillo) and GSE164690 (Kürten) from GEO.

---

## Pipeline

We run two pseudobulking strategies. `bulk_before` pseudobulks expression and then runs
scFEA. This is the primary arm and every reported number comes from it. `bulk_after` runs
scFEA per cell and then pseudobulks. This is the robustness check. We treat a finding as
robust only when both arms agree.

The analysis unit is **patient x tissue x cell subtype**, and needs at least 10 cells.

### Stage 0: imputation (cluster, per major cell type)

| Script | Produces |
|---|---|
| `cloud_cluster/all_sites/shared/00a_export_normalized_sparse.R` | `normalized.mtx`, `genes.txt`, `barcodes.txt`, `cell_metadata.csv` |
| `cloud_cluster/all_sites/shared/00b_magic_impute.py --solver exact` | `imputed_expression.csv` (661 metabolic genes x cells) |
| `cloud_cluster/all_sites/shared/00c_transpose_imputed.py` | scFEA orientation |

MAGIC runs once per major cell type. The graph is built on all genes, and only the 661
scFEA module genes are returned. Both arms read the same imputed matrices, so imputation
is identical and the arms differ only in aggregation order.

T is the one to watch. It needs a largemem node and a 6 hour job, so it runs through
`00b_magic_t.slurm` and the array in `00b_magic_impute.slurm` covers indices 2-9.

`submit_all.sh` queues Stages 0 and 1 in dependency order. Run it with no arguments for a
dry run.

### Stage 1: pseudobulk and flux (cluster)

| Script | Produces |
|---|---|
| `cloud_cluster/all_sites/bulk_before/01_export_bulk_before.R` | `train_expression.csv`, `metadata_units.csv` |
| `cloud_cluster/all_sites/bulk_before/02_scfea_bulk_before.slurm` | module flux per unit |
| `cloud_cluster/all_sites/bulk_after/01_export_bulk_after.R`, then `02_train`, then `03_predict` | per-cell flux (secondary arm) |

### Stage 2: assemble

Everything from here runs locally. `run_all.sh` does it in order.

```
01_build_combined_wide.R        ->  cw_finaltype.rds
01d_floor_and_mark_comparable.R ->  negative flux floored at zero
01e_fix_unit_duplicates.R       ->  canonical units (bulk_after only)
03_qc_flux.R                    ->  usable_modules.rds
01b_build_combined_wide_after.R ->  same, secondary arm
```

### Stage 3: analysis

```
04_specialization.R           ->  specialization_all.csv        (feeds 07b)
05_progression.R              ->  progression_all.csv           (Fig 2A)
05c_progression_meta.R        ->  meta_all.csv                  (Fig 2B, Tables S1-S3)
06_hpv.R                      ->  hpv_all.csv                   (Table S4)
07b_grounding_scalefree.R     ->  gene_level_scalefree.csv      (Fig 5A)
07c_grounding_final.R         ->  module_grounding_final.csv
07e_null_max_matched.R        ->  null_max_matched_thresholds.csv
07g_patch_matched_threshold.R ->  rewrites module_grounding_final.csv IN PLACE
make_progression_deg_edger.R  ->  edgeR differential expression
```

> **Order matters here.** `07e_null_max_matched.R` computes the matched empirical null
> (the 99.9th percentile of the maximum of *k* correlations, where *k* is the number of
> evaluable genes in the module). `07g_patch_matched_threshold.R` then writes those
> thresholds into `module_grounding_final.csv`, replacing the single-draw values from
> `07c_grounding_final.R`. **Both must run before any figure.** If you skip them the
> figures use the old threshold without telling you. The original file is kept as
> `module_grounding_final_singledraw_backup.csv`.

### Stage 4: external cohorts

```
30_tcga_signature.R        ->  cox_univariate.csv, per_patient_scores.csv   (Fig 6A/B, Fig S)
33_tcga_deconvolution.R    ->  cox_adjusted.csv                             (Fig 6C)
31_cptac_protein.py        ->  tumor_vs_normal_protein.csv,
                               rna_protein_concordance.csv,
                               mct4_protein_correlations.csv                (Fig 5C/D)
32a_build_cillo.R, 32b_build_kurten.R  ->  cohort objects
32_independent_cohort.R    ->  C2_isoform_specificity.csv                   (Fig 5B)
```

`31_cptac_protein.py` needs the `cptac` environment, not R, so `run_all.sh` does not call
it. Run it before Figure 5.

### Stage 5: figures and tables

| Output | Script |
|---|---|
| Figure 1 | `fig1_atlas_overview.R` |
| Figure 2 | `fig2_allocation_panels.R` |
| Figure 3 | `fig3_reorganisation_panels.R` |
| Figure 4 | `fig4_progression_all_panels.R` |
| Figure 5 | `fig5_transporter_all_panels.R` |
| Figure 6 | `fig6_clinical_panels.R` |
| Figure S (hypoxia) | `figS_hypoxia.R` |
| Tables S1, S2 | `make_supplementary_tables.R` |
| Table S3 | `table_S3_survivors.R` |
| Table S4 | `make_table_S4_hpv.R` |
| Module ranking | `rank_modules_by_effect.R` |

Every figure reads from ten files and nothing else: `cw_finaltype.rds`,
`usable_modules.rds`, `progression_all.csv`, `meta_all.csv`,
`module_grounding_final.csv`, `gene_level_scalefree.csv`, `train_expression.csv`,
`hpv_all.csv`, and the CPTAC and TCGA output sets.

---

## supporting_analyses/

These are not on the path to a figure, but each one produced a number or a decision we
report in the paper.

| Script | Supports |
|---|---|
| `02_gate_sweep.R` | the 10-cell gate. Same coverage as 8 with 350 fewer degenerate units, and 15 costs two progression comparisons |
| `01c_negative_flux_diagnostic.R` | negative flux in 40 of 168 modules (secondary) against 6 of 168 (primary) |
| `01g_dominance_vs_agreement.R` | no single myeloid subtype dominates the compartment shift (n = 4-8) |
| `07f_impute_vs_raw.R` | the flux-enzyme correlation is not an artifact of imputation |
| `09_compositional_vs_intrinsic.R`, `09d_compartment_composition_test.R` | reprogramming against changing cell composition |
| `34c_exhaustion_partial_correlation.R` | T cell dysfunction association |
| `check_*.R` | targeted checks quoted in the text |

---

## Modified scFEA

We ran a modified copy of scFEA against upstream commit
`4c1fb76d52f07bafad84ce7686ad7c3acfcf0126` (9 November 2022). The changes are in
`cloud_cluster/scfea_modifications.patch`:

- `--mode {train,predict}` and `--model_file`, so a model trained on a subset can be
  applied to the full atlas
- `.cpu()` before `.numpy()`, because PyTorch tensors on a GPU cannot convert directly
- output paths respect `--res_dir` instead of a hardcoded `./output/`

To rebuild it:

```bash
git clone https://github.com/changwn/scFEA.git
cd scFEA && git checkout 4c1fb76d52f07bafad84ce7686ad7c3acfcf0126
git apply --ignore-whitespace /path/to/scfea_modifications.patch
```

We left the hyperparameters at the code defaults (`LEARN_RATE` 0.008, `LAMB_BA` 1,
`LAMB_NG` 1, `LAMB_CELL` 1, `LAMB_MOD` 1e-2, `torch.manual_seed(16)`, `--train_epoch
100`). They are set explicitly in every submission script and are the same in both arms.

The MIT license covers the analysis code here. `scfea_modifications.patch` is a
modification of scFEA and carries scFEA's own academic, non-commercial license.

---

## Data

Source single-cell cohorts: HRA004648 (He), GSE181919 (Kim), GSE198315 (Koh),
GSE234933 (Pittet), GSE182227 (Tirosh). HRA004648 is controlled access and requires an
approved request.

External validation cohorts: GSE139324 (Cillo), GSE164690 (Kürten).

Proteomics: CPTAC-HNSCC through the `cptac` Python package. Survival and bulk expression:
TCGA-HNSC through cBioPortal.

---

## Environment

| File | What it covers |
|---|---|
| `envs/local_r_sessionInfo.txt` | R session for Stages 2 to 5 |
| `envs/local_cptac.yml` | `cptac_env`, Python 3.10, for `31_cptac_protein.py`. cptac 1.5.14, umich pipeline, 182 samples x 12,224 proteins |
| `envs/cluster_scfea.yml` | `scfea_37`, Python 3.7 with magic-impute and torch, for MAGIC and scFEA |
| `envs/cluster_r.yml` | `covarnet_r`, for the R steps on the cluster (`00a`, `01_export_bulk_before`, `01_export_bulk_after`) |

`envs/cluster_r.yml` is a snapshot of a shared lab environment, so it records what we ran
in rather than something we controlled.

A conda export records what was installed, not what the code needs. If rebuilding from
one hits a solver conflict, use the package list under Requirements instead.

## Known gaps

- `jt_test()` (Jonckheere-Terpstra) is written out separately in `05`, `05c` and `09`
  instead of being sourced from one place.
