#!/usr/bin/env python3
# ==============================================================================
# 31_cptac_protein_check.py  --  does the MCT4 finding hold at the PROTEIN level?
#
#   python 31_cptac_protein_check.py
#
# WHY THIS IS THE IHC SUBSTITUTE.
# Every result so far is RNA-derived: scFEA flux comes from expression, and the
# enzyme grounding correlates that flux against the same expression. Internally
# consistent, but not independent. CPTAC-HNSCC has ~110 tumors with MATCHED
# proteomics and RNA, so it can ask a question our data structurally cannot:
# is the MCT4 protein actually there, and does it track the niche markers?
#
# STEP 0 IS A GO/NO-GO. Membrane transporters are notoriously under-detected in
# shotgun proteomics -- they are hydrophobic, low-abundance, and tryptic peptides
# are scarce. SLC16A3 may simply be absent from the protein table. That is not a
# negative result about the biology, it is a limit of the assay, and we need to
# know it in one minute rather than after building an analysis on it.
#
# IF SLC16A3 IS PRESENT, we test:
#   1. protein-protein: does SLC16A3 track SPP1 / MMP11 protein across tumors?
#   2. RNA-protein concordance for SLC16A3 (is the transcript a fair proxy?)
#   3. tumor vs normal-adjacent: is MCT4 protein higher in tumor?
#   4. the same for the supporting genes (LDHA, HK2) as internal controls
#
# CAVEAT THAT SURVIVES ANY RESULT: CPTAC is BULK protein. A correlation between
# SLC16A3 and SPP1 protein is consistent with "MCT4 is on SPP1+ macrophages" but
# also with "tumors with more macrophages have more of both". Bulk cannot
# separate those. This corroborates; it does not prove the cell-intrinsic claim.
# ==============================================================================
import sys, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from scipy import stats

OUT = "/Users/Tommy/dechenli_lab_local/scfea26/results/cptac"
import os; os.makedirs(OUT, exist_ok=True)

TARGET   = "SLC16A3"                    # MCT4 -- the claim
NICHE    = ["SPP1", "MMP11", "CD68"]    # niche identity
SUPPORT  = ["LDHA", "HK2", "SLC16A1"]   # supporting / contrast (MCT1)
NEGCTRL  = ["ACTB", "B2M", "GAPDH"]     # should NOT be specifically informative
ALL_GENES = [TARGET] + NICHE + SUPPORT + NEGCTRL

# ------------------------------------------------------------------ load ----
def load_cptac():
    """cptac's API changed at v1.0; try the modern form, then the legacy one."""
    try:
        import cptac
    except ImportError:
        sys.exit("pip install cptac   (then re-run)")

    print(f"cptac version: {getattr(cptac, '__version__', 'unknown')}")
    # --- modern API (>=1.0) ---
    try:
        hn = cptac.Hnscc()
        srcs = hn.list_data_sources()
        print("data sources available:")
        print(srcs.to_string() if hasattr(srcs, "to_string") else srcs)
        prot = None
        for src in ("umich", "bcm", "broad", "washu"):
            try:
                prot = hn.get_proteomics(source=src); print(f"  proteomics from: {src}")
                break
            except Exception:
                continue
        if prot is None:
            raise RuntimeError("no proteomics source worked")
        try:
            rna = hn.get_transcriptomics(source="bcm")
        except Exception:
            rna = None
        return hn, prot, rna
    except Exception as e:
        print(f"modern API failed ({e}); trying legacy...")

    # --- legacy API (<1.0) ---
    try:
        cptac.download(dataset="hnscc")
    except Exception:
        pass
    hn = cptac.Hnscc()
    prot = hn.get_proteomics()
    try:    rna = hn.get_transcriptomics()
    except Exception: rna = None
    return hn, prot, rna

hn, prot, rna = load_cptac()

# columns are sometimes a MultiIndex (gene, database_id) -- flatten to gene
if isinstance(prot.columns, pd.MultiIndex):
    prot.columns = prot.columns.get_level_values(0)
if rna is not None and isinstance(rna.columns, pd.MultiIndex):
    rna.columns = rna.columns.get_level_values(0)

print(f"\nproteomics matrix: {prot.shape[0]} samples x {prot.shape[1]} proteins")

# ================= STEP 0: THE GO / NO-GO =====================================
print("\n" + "="*62)
print("STEP 0 -- IS MCT4 (SLC16A3) EVEN DETECTED IN THE PROTEOME?")
print("="*62)
present = {g: (g in prot.columns) for g in ALL_GENES}
for g, ok in present.items():
    if ok:
        v = prot[g]
        if isinstance(v, pd.DataFrame): v = v.iloc[:, 0]
        print(f"  [{'ok' if ok else '--'}] {g:<9} detected in "
              f"{v.notna().mean()*100:5.1f}% of samples")
    else:
        print(f"  [--] {g:<9} NOT in the protein table")

if not present[TARGET]:
    print(f"\n  {TARGET} IS NOT MEASURED in CPTAC-HNSCC proteomics.")
    print("  This is an ASSAY LIMIT, not evidence against the finding: membrane")
    print("  transporters are hydrophobic, low-abundance, and yield few tryptic")
    print("  peptides, so they drop out of shotgun proteomics routinely.")
    print("\n  -> Report as: 'MCT4 protein was not quantified in available HNSCC")
    print("     proteomic data, so protein-level validation requires targeted")
    print("     assays (IHC / PRM).' Then move to the independent scRNA cohort,")
    print("     which is the stronger validation anyway.")
    pd.DataFrame([{"gene": g, "detected": ok} for g, ok in present.items()]) \
      .to_csv(f"{OUT}/detection_check.csv", index=False)
    sys.exit(0)

print(f"\n  {TARGET} IS DETECTED -- proceeding to the correlation tests.")

# ------------------------------------------------------- helper --------------
def col(df, g):
    if g not in df.columns: return None
    v = df[g]
    return v.iloc[:, 0] if isinstance(v, pd.DataFrame) else v

def spearman(a, b):
    m = a.notna() & b.notna()
    if m.sum() < 20: return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(a[m], b[m])
    return r, p, int(m.sum())

# ============ 1. does MCT4 protein track the niche markers? ==================
print("\n" + "="*62)
print("1 -- MCT4 PROTEIN vs NICHE MARKER PROTEIN (across tumors)")
print("="*62)
tgt = col(prot, TARGET)
rows = []
for g in NICHE + SUPPORT + NEGCTRL:
    o = col(prot, g)
    if o is None: continue
    r, p, n = spearman(tgt, o)
    grp = ("niche" if g in NICHE else "support" if g in SUPPORT else "neg control")
    rows.append(dict(gene=g, group=grp, rho=r, p=p, n=n))
cor_tab = pd.DataFrame(rows).sort_values("rho", ascending=False)
print(cor_tab.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
cor_tab.to_csv(f"{OUT}/mct4_protein_correlations.csv", index=False)
print("\n  READ: niche markers (SPP1/MMP11/CD68) clearly above the negative")
print("  controls (ACTB/B2M) is the supportive pattern. If ACTB tracks MCT4 as")
print("  well as SPP1 does, the correlation is loading/abundance, not biology.")

# ============ 2. RNA-protein concordance ====================================
if rna is not None:
    print("\n" + "="*62)
    print("2 -- RNA vs PROTEIN concordance (is the transcript a fair proxy?)")
    print("="*62)
    shared = prot.index.intersection(rna.index)
    print(f"  samples with both: {len(shared)}")
    rows = []
    for g in [TARGET] + NICHE + SUPPORT:
        pv, rv = col(prot, g), col(rna, g)
        if pv is None or rv is None: continue
        r, p, n = spearman(pv.loc[shared], rv.loc[shared])
        rows.append(dict(gene=g, rho_rna_protein=r, p=p, n=n))
    conc = pd.DataFrame(rows)
    print(conc.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    conc.to_csv(f"{OUT}/rna_protein_concordance.csv", index=False)
    print("\n  READ: a decent positive rho for SLC16A3 means our RNA-based flux")
    print("  is tracking something that exists as protein -- the key assumption")
    print("  behind every scFEA result. A rho near zero would be a real caveat.")

# ============ 3. tumor vs normal-adjacent ===================================
print("\n" + "="*62)
print("3 -- TUMOR vs NORMAL-ADJACENT (is MCT4 protein higher in tumor?)")
print("="*62)
try:
    clin = hn.get_clinical(source="mssm")
except Exception:
    try: clin = hn.get_clinical()
    except Exception: clin = None

tumor_mask = None
if clin is not None:
    for c in ("Sample_Tumor_Normal", "tumor_normal", "Type", "sample_type"):
        if c in clin.columns:
            s = clin[c].reindex(prot.index).astype(str).str.lower()
            tumor_mask = s.str.contains("tumor")
            break
if tumor_mask is None:                       # CPTAC normals often end in ".N"
    tumor_mask = ~prot.index.astype(str).str.endswith(".N")

print(f"  tumor: {int(tumor_mask.sum())} | normal: {int((~tumor_mask).sum())}")
if (~tumor_mask).sum() >= 10:
    rows = []
    for g in [TARGET] + NICHE + SUPPORT:
        v = col(prot, g)
        if v is None: continue
        t, nrm = v[tumor_mask].dropna(), v[~tumor_mask].dropna()
        if len(t) < 10 or len(nrm) < 10: continue
        u, p = stats.mannwhitneyu(t, nrm)
        rows.append(dict(gene=g, tumor_med=t.median(), normal_med=nrm.median(),
                         diff=t.median()-nrm.median(), p=p))
    tn = pd.DataFrame(rows).sort_values("diff", ascending=False)
    print(tn.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    tn.to_csv(f"{OUT}/tumor_vs_normal_protein.csv", index=False)
else:
    print("  too few normal-adjacent samples for this comparison")

print(f"\n-> {OUT}")
print("\nCAVEAT for the writeup: CPTAC is BULK protein. MCT4 tracking SPP1 is")
print("consistent with MCT4 being ON SPP1+ macrophages, and equally consistent")
print("with macrophage-rich tumors simply having more of both. Bulk cannot")
print("separate them -- this corroborates the finding, it does not prove the")
print("cell-intrinsic part.")
