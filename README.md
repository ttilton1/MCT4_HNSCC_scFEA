# MCT4-associated lactate export in HNSCC myeloid cells

Analysis code for:

> Tilton T, Park YM, Zhu H, Zhang C, Lin D-C. *MCT4-Associated Lactate Export Defines
> Myeloid Metabolic Reprogramming during Head and Neck Cancer Progression.*
> bioRxiv 2026. https://doi.org/10.64898/2026.10.13.758108

Metabolic flux is inferred with [scFEA](https://github.com/changwn/scFEA) across 168
modules (M1–M168) of the M171 network, in 645,351 cells from 114 patients across five
published HNSCC cohorts. The analysis is organised around **metabolic allocation**, defines as a 
module's flux as a fraction of its unit's total flux so we can compare across cell types. Using
raw fluxes would not let us compare across major cell types that each have their own
flux model to fit to.

---

## Setup

Scripts resolve the data directory from an environment variable:

```bash
export SCFEA_ROOT=/path/to/your/scfea_data
```

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

**Requirements.** R ≥ 4.2 with dplyr, tidyr, readr, data.table, Matrix, ggplot2,
patchwork, ggtext, survival, edgeR, Seurat. Python 3.7 with `magic-impute` and torch for
the cluster steps, and a separate environment with `cptac==1.5.14` for the proteomics.
`environment.yml` and `sessionInfo.txt` record the versions actually used.

**External data not included here.** The five source cohorts (He, Kim, Koh, Pittet,
Tirosh) as a harmonised Seurat object; TCGA-HNSC from cBioPortal; CPTAC-HNSCC through
the `cptac` package; GSE139324 (Cillo) and GSE164690 (Kürten) from GEO.

---

## Pipeline

Two pseudobulking strategies run in parallel. `bulk_before` — pseudobulk expression,
then scFEA — is the **primary** arm and the source of every reported number.
`bulk_after` — scFEA per cell, then pseudobulk — is the robustness check. A finding is
treated as robust only when both agree.

The analysis unit is **patient × tissue × cell subtype**, requiring ≥10 cells.

### Stage 0 — imputation (cluster, per major cell type)

| Script | Produces |
|---|---|
| `cloud_cluster/all_sites/shared/00a_export_normalized_sparse.R` | `normalized.mtx`, `genes.txt`, `barcodes.txt`, `cell_metadata.csv` |
| `cloud_cluster/all_sites/shared/00b_magic_impute.py --solver exact` | `imputed_expression.csv` (661 metabolic genes × cells) |
| `cloud_cluster/all_sites/shared/00c_transpose_imputed.py` | scFEA orientation |

MAGIC runs once per major cell type, with the neighbourhood graph built on all genes and
only the 661 scFEA module genes returned. Both arms read the same imputed matrices, so
imputation is held identical and the arms differ only in aggregation order. The T cell
type needs a largemem node — see `00b_magic_t.slurm`.

### Stage 1 — pseudobulk and flux (cluster)

| Script | Produces |
|---|---|
| `cloud_cluster/all_sites/bulk_before/01_export_bulk_before.R` | `train_expression.csv`, `metadata_units.csv` |
| `cloud_cluster/all_sites/bulk_before/02_scfea_bulk_before.slurm` | module flux per unit |
| `cloud_cluster/all_sites/bulk_after/01_export_bulk_after.R` → `02_train` → `03_predict` | per-cell flux (secondary arm) |

### Stage 2 — assemble

```
01_build_combined_wide.R        →  cw_finaltype.rds
01d_floor_and_mark_comparable.R →  negative flux floored at zero
01e_fix_unit_duplicates.R       →  canonical units (bulk_after only)
03_qc_flux.R                    →  usable_modules.rds
01b_build_combined_wide_after.R →  same, secondary arm
```

### Stage 3 — analysis

```
04_specialization.R       →  specialization_all.csv        (feeds 07b)
05_progression.R          →  progression_all.csv           (Fig 2A)
05c_progression_meta.R    →  meta_all.csv                  (Fig 2B, Tables S1–S3)
06_hpv.R                  →  hpv_all.csv                   (Table S4)
07b_grounding_scalefree.R →  gene_level_scalefree.csv      (Fig 5A)
07c_grounding_final.R     →  module_grounding_final.csv
07e_null_max_matched.R    →  null_max_matched_thresholds.csv
07c_patch_matched.R       →  rewrites module_grounding_final.csv IN PLACE
make_progression_deg_edger.R → edgeR differential expression
```

> **Ordering matters here.** `07e_null_max_matched.R` computes the matched empirical
> null — the 99.9th percentile of the maximum of *k* correlations, where *k* is the
> number of evaluable genes in the module — and `07c_patch_matched.R` writes those
> thresholds into `module_grounding_final.csv`, replacing the single-draw values from
> `07c_grounding_final.R`. **Both must run before any figure.** Skip them and the
> figures silently use the superseded threshold. The original is preserved as
> `module_grounding_final_singledraw_backup.csv`.

### Stage 4 — external cohorts

```
30_tcga_signature.R        →  cox_univariate.csv, per_patient_scores.csv   (Fig 6A/B, Fig S)
33_tcga_deconvolution.R    →  cox_adjusted.csv                             (Fig 6C)
31_cptac_protein_check.py  →  tumor_vs_normal_protein.csv, rna_protein_concordance.csv,
                              mct4_protein_correlations.csv                (Fig 5C/D)
32a_build_cillo.R, 32b_build_kurten.R  →  cohort objects
32_independent_cohort.R    →  C2_isoform_specificity.csv                   (Fig 5B)
```

### Stage 5 — figures and tables

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

Every figure reads from ten intermediate files and nothing else:
`cw_finaltype.rds`, `usable_modules.rds`, `progression_all.csv`, `meta_all.csv`,
`module_grounding_final.csv`, `gene_level_scalefree.csv`, `train_expression.csv`,
`hpv_all.csv`, and the CPTAC and TCGA output sets.

---

## supporting_analyses/

Not on the path to a figure, but each produced a number or a decision reported in the
paper:

| Script | Supports |
|---|---|
| `02_gate_sweep.R` | the 10-cell gate: identical coverage to 8, 350 fewer degenerate units; 15 costs two progression comparisons |
| `01c_negative_flux_diagnostic.R` | negative flux in 40 of 168 modules (secondary) vs 6 of 168 (primary) |
| `01g_dominance_vs_agreement.R` | no single myeloid subtype dominates the compartment shift (n = 4–8) |
| `07f_impute_vs_raw.R` | the flux–enzyme correlation is not an artefact of imputation |
| `09_compositional_vs_intrinsic.R`, `09d_compartment_composition_test.R` | reprogramming vs changing cell composition |
| `34c_exhaustion_partial_correlation.R` | T-cell dysfunction association |
| `check_*.R` | targeted checks quoted in the text |

---

## Modified scFEA

Flux was estimated with a modified copy of scFEA, against upstream commit
`4c1fb76d52f07bafad84ce7686ad7c3acfcf0126` (9 November 2022). The changes are in
`cloud_cluster/scfea_modifications.patch`:

- `--mode {train,predict}` and `--model_file`, so a model trained on a subset can be
  applied to the full atlas
- `.cpu()` before `.numpy()`, since PyTorch tensors on a GPU cannot convert directly
- output paths respect `--res_dir` rather than a hardcoded `./output/`

Apply with:

```bash
git clone https://github.com/changwn/scFEA.git
cd scFEA && git checkout 4c1fb76d52f07bafad84ce7686ad7c3acfcf0126
git apply --ignore-whitespace /path/to/scfea_modifications.patch
```

Hyperparameters were left at the code defaults (`LEARN_RATE` 0.008, `LAMB_BA` 1,
`LAMB_NG` 1, `LAMB_CELL` 1, `LAMB_MOD` 1e-2, `torch.manual_seed(16)`,
`--train_epoch 100`), set explicitly in every submission script and identical across
both arms.

---

## Known gaps

- `jt_test()` (Jonckheere–Terpstra) is defined separately in `05`, `05c` and `09` rather
  than sourced from one place.

### Environment

`environment.yml` and `sessionInfo.txt` record the R environment. Two Python
environments are not captured there:

- **Proteomics** — `cptac==1.5.14`, using the umich pipeline (182 samples × 12,224
  proteins). Required by `31_cptac_protein.py`.
- **Cluster** — `scfea_37` (Python 3.7) with `magic-impute` and torch, for Stages 0–1.
