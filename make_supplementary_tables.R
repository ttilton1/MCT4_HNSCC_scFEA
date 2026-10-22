#!/usr/bin/env Rscript
# ==============================================================================
# make_supplementary_tables.R  --  Supplementary Tables S2 and S3 as TSV
#
#   Rscript make_supplementary_tables.R --arm bulk_before
#
# S2  all 19 progression survivors, with per-cohort effects and criterion status.
# S3  paired pseudobulk differential expression for myeloid, fibroblast and T
#     cells, restricted to genes that moved consistently, since the full tables
#     run to 9,000 rows and belong in the repository rather than a supplement.
#
# Both are written tab-separated so they paste straight into a spreadsheet.
#
# Requires 05c_progression_meta.R, 07c_grounding_final.R and
# make_progression_deg.R to have run.
#
# Output -> results/all_sites/<ARM>/supplementary/
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

OUT <- file.path(PATHS$results, "supplementary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ------------------------------- TABLE S2 -----------------------------------
meta <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
                 show_col_types = FALSE)
gr   <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
                 show_col_types = FALSE) %>%
        select(celltype, module, best_gene, best_rho, null_bar, GROUNDED)

s2 <- meta %>% filter(survivor) %>%
  left_join(gr, by = c("celltype","module")) %>%
  mutate(criterion = case_when(is.na(GROUNDED) ~ "not tested",
                               GROUNDED        ~ "meets",
                               TRUE            ~ "does not meet"),
         in_paper = ifelse(unit == "DC_C2_cDC2" & module %in% c("M_78","M_6"), "yes", "")) %>%
  arrange(factor(criterion, levels = c("meets","does not meet","not tested")),
          desc(abs(delta_meta))) %>%
  transmute(`cell type` = celltype, subtype = unit, module, reaction,
            `change (meta)` = round(delta_meta, 4),
            `change He` = round(delta_He, 4), `change Kim` = round(delta_Kim, 4),
            `patients He` = n_pt_He, `patients Kim` = n_pt_Kim,
            FDR = signif(fdr, 3),
            `best catalytic gene` = best_gene, `rho` = round(best_rho, 3),
            `null threshold` = round(null_bar, 3),
            `flux-enzyme criterion` = criterion,
            `carried forward` = in_paper)

write_tsv(s2, file.path(OUT, "TABLE_S2_progression_survivors.tsv"))
say("=== Table S2: %d progression survivors ===", nrow(s2))
say("  meeting the criterion: %d | not: %d | not tested: %d",
    sum(s2$`flux-enzyme criterion` == "meets"),
    sum(s2$`flux-enzyme criterion` == "does not meet"),
    sum(s2$`flux-enzyme criterion` == "not tested"))

# ------------------------------- TABLE S3 -----------------------------------
DEGDIR <- file.path(ROOT, "results", "deg")
lins <- c(Myeloid = "myeloid", Fibroblast = "fibroblast", `T cells` = "t")
s3 <- lapply(names(lins), function(nm) {
  f <- file.path(DEGDIR, sprintf("DEG_%s_normal_vs_carcinoma_edger.csv", lins[[nm]]))
  if (!file.exists(f)) { say("  missing: %s", f); return(NULL) }
  d <- read_csv(f, show_col_types = FALSE)
  d %>% filter(n_rising == n_patients | n_rising == 0, FDR < 0.10) %>%
    arrange(desc(abs(logFC))) %>%
    transmute(`cell type` = nm, gene,
              log2FC = round(logFC, 3),
              logCPM = round(logCPM, 2),
              `patients rising` = sprintf("%d/%d", n_rising, n_patients),
              `p` = signif(PValue, 3), FDR = signif(FDR, 3))
}) %>% bind_rows()

write_tsv(s3, file.path(OUT, "TABLE_S3_differential_expression.tsv"))
say("\n=== Table S3: %d genes moving consistently at FDR < 0.10 ===", nrow(s3))
print(as.data.frame(s3 %>% count(`cell type`, name = "genes")), row.names = FALSE)
say("\n  restricted to genes changing in every paired patient (or none) at FDR < 0.10.")
say("  Full edgeR tables (12,028 to 16,947 genes) belong in the repository.")
say("\n-> %s", OUT)
