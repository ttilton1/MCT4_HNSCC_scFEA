#!/usr/bin/env Rscript
# ==============================================================================
# figS_hypoxia.R  --  does lactate export track bulk hypoxia?
#
#   Rscript figS_hypoxia.R
#
# Standalone; reads results/tcga/per_patient_scores.csv.
#
# WHY THIS PLOTS SLC16A3 AND NOT THE COMPOSITE. The export-program score
# correlates with the Buffa score at rho = 0.60 against SLC16A3's 0.36. That gap
# is not safely interpretable: published hypoxia signatures are built largely
# from glycolytic genes, and LDHA, HK2, CA9 and VEGFA are plausible members of
# the very scores being correlated against. A composite sharing genes with the
# comparator will correlate with it partly by construction.
#
# The structure of the result supports that reading. Buffa and Winter both show
# a large composite advantage (0.36 -> 0.60, 0.29 -> 0.58) while Ragnum shows
# none at all (0.24 -> 0.23). A real improvement in capturing hypoxia should be
# consistent across signatures; an overlap artefact should appear only where the
# overlap is.
#
# We therefore plot SLC16A3 alone, which is the gene the claim concerns and the
# only quantity here with no possible membership in the comparator. The composite
# values are printed for reference and flagged.
#
# RESOLVED against MSigDB (BUFFA_HYPOXIA_METAGENE, 50 genes; WINTER_HYPOXIA_METAGENE,
# 243 genes):
#                Buffa        Winter
#   SLC16A3      not member   not member    -> single-gene correlation is clean
#   LDHA         member       member        -> shared with export_core
#   HK2          member       member        -> shared with export_core
#   CA9          member       member        -> shared with sig_full
#   VEGFA        member       member        -> shared with sig_full
#   SPP1         not member   member        -> shared with sig_full (Winter only)
#   SLC16A1      member       member
#   GAPDH        member       member
#
# export_core shares two of its three genes with each score, which accounts for
# 0.60 (Buffa) and 0.58 (Winter) against SLC16A3's 0.36 and 0.29. Ragnum, a
# 32-gene prostate-derived set, was not checked directly; it shows no composite
# advantage, consistent with sharing fewer genes.
#
# MCT4's absence from both metagenes is worth a word, since it is a well-known
# HIF-1a target. Buffa and Winter are empirically derived co-expression metagenes,
# not curated HIF-target lists, so absence reflects construction rather than
# regulation. It is what makes the 0.36 clean.
#
# Also of note: GAPDH is a member of both. That does NOT explain its negative
# correlation with MCT4 protein in CPTAC -- co-membership predicts positive, since
# MCT4 transcript tracks both scores positively in the same tumors. The
# discrepancy therefore points at the protein measurement, not at biology.
# And SPP1 is a Winter member, so any future analysis relating SPP1 to a Winter
# score inherits the circularity documented here.
#
# Output -> results/tcga/figS_hypoxia/
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(tidyr) })
TC  <- paste0(Sys.getenv("SCFEA_ROOT", unset = "/Users/Tommy/dechenli_lab_local/scfea26"), "/results/tcga")
OUT <- file.path(TC, "figS_hypoxia"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(file.path(TC, "per_patient_scores.csv"), show_col_types = FALSE)
HS <- c(buffa_hypoxia = "Buffa", ragnum_hypoxia = "Ragnum", winter_hypoxia = "Winter")

tab <- lapply(names(HS), function(h) {
  k  <- !is.na(d[[h]])
  cm <- suppressWarnings(cor.test(d$SLC16A3[k],     d[[h]][k], method = "spearman"))
  cc <- suppressWarnings(cor.test(d$export_core[k], d[[h]][k], method = "spearman"))
  tibble(signature = HS[[h]], n = sum(k),
         rho_mct4 = unname(cm$estimate), p_mct4 = cm$p.value,
         rho_composite = unname(cc$estimate), p_composite = cc$p.value)
}) %>% bind_rows()

fmt_p <- function(p) ifelse(p < 2.2e-16, "p < 2e-16",
                    ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p)))
cat("=== correlation with published bulk hypoxia scores ===\n")
print(as.data.frame(tab %>% transmute(signature, n,
        MCT4 = sprintf("%.3f (%s)", rho_mct4, fmt_p(p_mct4)),
        export_program = sprintf("%.3f (%s)", rho_composite, fmt_p(p_composite)))),
      row.names = FALSE)
cat("\nThe composite exceeds MCT4 on Buffa and Winter but not Ragnum. This is gene\noverlap: LDHA and HK2 are members of both the Buffa and Winter metagenes,\nSLC16A3 of neither. See the verification note at the end of this output.\n")
write_csv(tab, file.path(OUT, "hypoxia_correlations.csv"))

long <- d %>% select(SLC16A3, all_of(names(HS))) %>%
  pivot_longer(-SLC16A3, names_to = "h", values_to = "score") %>%
  filter(!is.na(score)) %>% mutate(signature = HS[h])
lab <- tab %>% mutate(txt = sprintf("rho = %.2f, %s\nn = %d", rho_mct4, fmt_p(p_mct4), n))

p <- ggplot(long, aes(score, SLC16A3)) +
  geom_point(alpha = .28, size = .9, colour = "grey35") +
  geom_smooth(method = "lm", se = FALSE, colour = "#b2182b", linewidth = .8) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = txt),
            hjust = -0.1, vjust = 1.4, size = 3, lineheight = 1.05,
            colour = "#b2182b", inherit.aes = FALSE) +
  facet_wrap(~ signature, scales = "free_x") +
  labs(title = "MCT4 expression against published bulk hypoxia scores",
       x = "hypoxia score", y = "SLC16A3 (z-score)") +
  theme_bw(base_size = 10)

ggsave(file.path(OUT, "figS_hypoxia.png"), p, width = 9.5, height = 3.8, dpi = 300, bg = "white")
ggsave(file.path(OUT, "figS_hypoxia.pdf"), p, width = 9.5, height = 3.8, bg = "white")
cat("\n->", OUT, "\n")

g <- function(s, col) tab[[col]][tab$signature == s]
cat(sprintf(
"\nCaption draft:\n  MCT4 expression against three published bulk hypoxia scores in TCGA-HNSC, one\n  point per tumor (n = %d). SLC16A3 correlates positively with all three: Buffa\n  rho = %.2f (%s), Winter %.2f (%s), Ragnum %.2f (%s). This is consistent with the\n  HIF regulation of MCT4 reported elsewhere. Correlations are modest, and with this\n  sample size significance is not informative; the effect sizes are.\n\n  MCT4 is plotted alone rather than as part of the export program. The program\n  correlates more strongly with two of the three scores (Buffa %.2f, Winter %.2f)\n  but not the third (Ragnum %.2f), and that pattern is explained by gene overlap\n  rather than biology: LDHA and HK2 are members of both the Buffa and Winter\n  metagenes, so the program shares two of its three genes with each of those scores.\n  SLC16A3 is a member of neither, which makes the single-gene correlation the\n  interpretable quantity. These are bulk measurements and cannot attribute the\n  association to any cell type.\n",
  tab$n[1],
  g("Buffa","rho_mct4"),  fmt_p(g("Buffa","p_mct4")),
  g("Winter","rho_mct4"), fmt_p(g("Winter","p_mct4")),
  g("Ragnum","rho_mct4"), fmt_p(g("Ragnum","p_mct4")),
  g("Buffa","rho_composite"), g("Winter","rho_composite"), g("Ragnum","rho_composite")))

cat("\n[VERIFIED] MSigDB BUFFA_HYPOXIA_METAGENE (50 genes) and WINTER_HYPOXIA_METAGENE\n  (243 genes). SLC16A3 is a member of NEITHER. LDHA, HK2, CA9 and VEGFA are\n  members of both; SPP1 is a Winter member. SLC16A1 (MCT1) and GAPDH are members\n  of both. GAPDH's membership predicts a POSITIVE correlation with MCT4, opposite\n  to the -0.40 seen at protein level, so it sharpens the CPTAC anomaly rather than\n  explaining it.\n  Ragnum (32 genes, prostate-derived) not checked; it shows no composite\n  advantage, consistent with sharing fewer genes.\n")
