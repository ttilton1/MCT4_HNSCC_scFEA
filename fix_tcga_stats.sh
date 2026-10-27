#!/usr/bin/env bash
# fix_tcga_stats.sh -- two statistical corrections to 30_tcga_signature.R.
#
#   cd /Users/Tommy/dechenli_lab_local/MCT4_HNSCC_scFEA
#   bash fix_tcga_stats.sh
#   git diff
#
# A. ICD-10 SITE MAPPING. case_when evaluates in order, and C01 (base of tongue)
#    matches the C00-C06 oral-cavity rule before it can reach the oropharynx
#    rule, so every base-of-tongue tumour was classified as oral cavity. Base of
#    tongue is oropharynx, and with tonsil is where HPV-positive disease
#    concentrates, so this shrinks and skews the within-oropharynx HPV test that
#    the manuscript reports as P = 0.94.
#
# B. PROPORTIONAL HAZARDS. The Cox models had no cox.zph() check. Adds one on
#    the reported model and writes it to cox_zph.csv.
set -o pipefail

cd "$(git rev-parse --show-toplevel 2>/dev/null)" || { echo "not in the repo"; exit 1; }
F=30_tcga_signature.R
[ -f "$F" ] || { echo "$F not found -- are you in MCT4_HNSCC_scFEA?"; exit 1; }

if [ -n "$(git status --porcelain $F)" ]; then
  echo "$F has uncommitted changes -- commit or stash first."; exit 1
fi

# ---- A. site mapping ---------------------------------------------------------
perl -0777 -pi -e 's{
  dh\s*<-\s*dh\s*%>%\s*mutate\(site\s*=\s*case_when\(.*?TRUE\s*~\s*"other"\)\)
}{dh <- dh %>% mutate(site = case_when(
      # C01 (base of tongue) is oropharynx, not oral cavity, and must be matched
      # BEFORE the C00-C06 rule or case_when captures it as oral cavity. With
      # tonsil (C09) it is where HPV-positive disease concentrates, so the order
      # determines who is in the within-oropharynx comparison below.
      grepl("C01|C09|C10", icd10) ~ "oropharynx",
      grepl("C0[0-6]", icd10)     ~ "oral_cavity",
      grepl("C32", icd10)         ~ "larynx",
      grepl("C12|C13", icd10)     ~ "hypopharynx",
      TRUE ~ "other"))}sx' "$F"

grep -q 'grepl("C01|C09|C10", icd10) ~ "oropharynx"' "$F" \
  && echo "ok   A. site mapping reordered" \
  || { echo "FAIL A. pattern did not match -- check the case_when block by hand"; exit 1; }

# ---- B. proportional hazards -------------------------------------------------
perl -0777 -pi -e 's{
(write_csv\(cox_tab,\s*file\.path\(OUT,\s*"cox_univariate\.csv"\)\)\n)
}{$1
# ---- proportional hazards check on the reported model ------------------------
# A Cox HR assumes the hazard ratio is constant over follow-up. cox.zph tests
# that; a small p for a term means its effect changes with time and the single
# HR is a summary of something non-constant.
.ph_d <- df %>% filter(!is.na(sig_full), !is.na(age))
.ph_f <- if (stage_ok) {
  .ph_d <- .ph_d %>% filter(!is.na(stage_simple))
  Surv(OS.time, OS) ~ sig_full + age + strata(stage_simple)
} else Surv(OS.time, OS) ~ sig_full + age
.ph <- cox.zph(coxph(.ph_f, data = .ph_d))
say("\\n=== proportional hazards (cox.zph), n = %d ===", nrow(.ph_d))
print(round(.ph$table, 4))
say("  a small p means that term\x27s effect changes over follow-up")
write_csv(tibble::rownames_to_column(as.data.frame(.ph$table), "term"),
          file.path(OUT, "cox_zph.csv"))
}sx' "$F"

grep -q "cox.zph" "$F" \
  && echo "ok   B. proportional hazards check added" \
  || echo "FAIL B. pattern did not match -- add cox.zph by hand"

echo
Rscript -e 'invisible(parse(commandArgs(TRUE)[1]))' "$F" \
  && echo "parses cleanly" || { echo "PARSE ERROR -- git checkout -- $F"; exit 1; }

cat <<'EOF'

Next:
  git diff
  Rscript 30_tcga_signature.R 2>&1 | tee /tmp/tcga_rerun.log

Then read these out of the log:
  grep -A3 "ICD-10 sub-sites"   /tmp/tcga_rerun.log
  grep -A2 "WITHIN OROPHARYNX"  /tmp/tcga_rerun.log
  grep -A6 "proportional"       /tmp/tcga_rerun.log

The oropharynx n should RISE. The manuscript reports P = 0.94 for that test, so
compare the new P against it before deciding whether the sentence needs changing.
EOF
