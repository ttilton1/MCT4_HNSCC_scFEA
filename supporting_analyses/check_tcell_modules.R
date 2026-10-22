#!/usr/bin/env Rscript
# ==============================================================================
# check_tcell_modules.R  --  what are T cells reallocating toward?
#
#   Rscript check_tcell_modules.R --arm bulk_before
#   Rscript check_tcell_modules.R --arm bulk_after
#
# WHY. Testing lactate export across lineages showed myeloid is the only one
# rising, and T cells are the only one falling significantly: -0.38 percentage
# points, 0 of 11 patients, p = 0.001, and -0.25, 0 of 11, p = 0.001 in the other
# arm. Shares sum to one, so if export falls something else rises. This asks
# what, applying the same tests the myeloid analysis used.
#
# WHAT TO WATCH FOR. Ten to eleven paired patients puts the Wilcoxon at a floor
# of p = 0.002, which most modules reach if the direction is consistent. Effect
# size among modules moving in every patient is the informative quantity, and
# only modules meeting the flux-enzyme criterion in T cells are interpretable.
# The T-cell threshold is 0.731, the highest of any cell type, so few will.
#
# A CAUTION SPECIFIC TO T CELLS. M_78 fails the criterion in T cells: the best
# correlate is SLC16A1 at 0.480, well under the threshold. So the export decline
# reported above is itself not interpretable by the paper's own standard. It is
# reported as a contrast to the myeloid rise, not as a claim about T-cell lactate
# handling, and anything found here inherits the same limitation unless the
# module clears the bar.
#
# Writes CSVs for upload alongside the T-cell DEG.
# ==============================================================================
source("config_all_sites.R")
suppressPackageStartupMessages({ library(dplyr) })

OUT <- file.path(PATHS$results, "tcell_modules")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
NL <- "NL"; CA <- c("CA","CA-E","CA-A","CA-c","CA-p")
u  <- readRDS(file.path(PATHS$results, "03_qc", "usable_modules.rds"))
mt <- read_csv(file.path(PATHS$results, "05c_progression_meta/finaltype/meta_all.csv"),
               show_col_types = FALSE) %>% distinct(module, reaction, supermodule)
gr <- read_csv(file.path(PATHS$results, "07c_grounding_final/module_grounding_final.csv"),
               show_col_types = FALSE) %>% filter(celltype == "T") %>%
      select(module, best_gene, best_rho, null_bar, GROUNDED)

cw <- as.data.frame(readRDS(file.path(PATHS$cw_dir, "cw_finaltype.rds"))) %>%
  filter(annot_4th == "T", n_cells >= GATES$min_cells,
         study %in% PROG_COHORTS, tissue %in% c(NL, CA)) %>%
  mutate(stage = ifelse(tissue %in% NL, "normal", "carcinoma"))
mc <- intersect(mcols_of(cw), u)
S  <- as.matrix(cw[, mc]); S <- S / rowSums(S)
say("=== ARM %s | T cells | %d units, %d modules ===", ARM, nrow(cw), length(mc))

res <- lapply(mc, function(m) {
  d <- cw %>% mutate(v = 100 * S[, m]) %>%
    group_by(study, patient, stage) %>%
    summarise(v = weighted.mean(v, n_cells), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = stage, values_from = v) %>%
    filter(!is.na(normal), !is.na(carcinoma))
  if (nrow(d) < 4) return(NULL)
  per <- sapply(split(d, d$study), function(z)
    if (nrow(z) >= 2) median(z$carcinoma) - median(z$normal) else NA)
  tibble(module = m, n = nrow(d),
         normal = median(d$normal), carcinoma = median(d$carcinoma),
         change = median(d$carcinoma) - median(d$normal),
         n_rising = sum(d$carcinoma > d$normal),
         He = unname(per["He"]), Kim = unname(per["Kim"]),
         p = suppressWarnings(wilcox.test(d$carcinoma, d$normal, paired = TRUE))$p.value)
}) %>% bind_rows() %>%
  mutate(FDR = p.adjust(p, "BH")) %>%
  left_join(mt, by = "module") %>% left_join(gr, by = "module") %>%
  mutate(criterion = ifelse(is.na(GROUNDED), "not tested",
                     ifelse(GROUNDED, "meets", "fails")),
         consistent = n_rising == n | n_rising == 0)

say("  at the Wilcoxon floor p = %.4f: %d of %d | meeting the criterion: %d",
    min(res$p), sum(res$p == min(res$p)), nrow(res), sum(res$criterion == "meets"))

say("\n=== modules RISING in every patient, by effect size ===")
up <- res %>% filter(n_rising == n, change > 0) %>% arrange(desc(change))
print(as.data.frame(up %>% slice_head(n = 12) %>%
  transmute(module, reaction = substr(reaction, 1, 32), change = round(change, 3),
            He = round(He, 3), Kim = round(Kim, 3),
            gene = best_gene, rho = round(best_rho, 3), criterion)), row.names = FALSE)
say("  %d modules rise in all patients | %d of those meet the criterion",
    nrow(up), sum(up$criterion == "meets"))

say("\n=== modules FALLING in every patient, by effect size ===")
dn <- res %>% filter(n_rising == 0) %>% arrange(change)
print(as.data.frame(dn %>% slice_head(n = 12) %>%
  transmute(module, reaction = substr(reaction, 1, 32), change = round(change, 3),
            He = round(He, 3), Kim = round(Kim, 3),
            gene = best_gene, rho = round(best_rho, 3), criterion)), row.names = FALSE)
say("  %d modules fall in all patients | %d of those meet the criterion",
    nrow(dn), sum(dn$criterion == "meets"))

say("\n=== the interpretable set: criterion-meeting and consistent ===")
ok <- res %>% filter(criterion == "meets", consistent) %>% arrange(desc(abs(change)))
if (nrow(ok)) {
  print(as.data.frame(ok %>% transmute(module, reaction = substr(reaction, 1, 34),
    normal = round(normal, 3), carcinoma = round(carcinoma, 3),
    change = round(change, 3), dir = ifelse(change > 0, "up", "down"),
    He = round(He, 3), Kim = round(Kim, 3), gene = best_gene,
    rho = round(best_rho, 3))), row.names = FALSE)
  say("  cohorts agreeing in sign: %d of %d", sum(sign(ok$He) == sign(ok$Kim), na.rm = TRUE), nrow(ok))
} else say("  none")

say("\n=== M_78 in T cells, for reference ===")
print(as.data.frame(res %>% filter(module == "M_78") %>%
  transmute(module, normal = round(normal,2), carcinoma = round(carcinoma,2),
            change = round(change,2), rising = sprintf("%d/%d", n_rising, n),
            p = signif(p,2), gene = best_gene, rho = round(best_rho,3), criterion)),
  row.names = FALSE)

write_csv(res %>% select(module, reaction, supermodule, n, normal, carcinoma,
                         change, n_rising, He, Kim, p, FDR,
                         best_gene, best_rho, criterion_threshold = null_bar, criterion),
          file.path(OUT, sprintf("tcell_flux_%s.csv", ARM)))
say("\n-> %s/tcell_flux_%s.csv", OUT, ARM)
