# Preregistration (written 2026-11-04, before any replicate was run)

This file fixes the predictions and the scoring rules **before** looking at results. If a result
contradicts a prediction, the README reports it as such; the predictions are not edited afterwards
(check `git log -- PREREGISTRATION.md`).

## Question
When a QTL acts through **allele dosage** in an autotetraploid, which genotype encodings can see it,
and what happens as the QTL sits in sequence that is more and more similar to a paralogous copy?

## Simulation (one replicate)
- Autotetraploid, 400 individuals; each homolog is a mosaic of 60 founder haplotypes (mean block 25 kb).
- chr1 2 Mb + chr2 0.5 Mb, random sequence (GC 36 %), SNP density 1/150 bp, substitutions only.
- Five traits, each = one major QTL (10 % of variance) + shared chr2 polygenic background (20 %, 100 loci) + noise (70 %).

| trait | QTL model | founder alt freq | where the QTL sits |
|---|---|---|---|
| U | additive dosage | 0.5 | unique sequence |
| CTRL | dominant (presence of ≥1 copy) | 0.1 | unique sequence |
| P95 | additive dosage | 0.5 | 20 kb segment with a paralogous copy at 95 % identity |
| P98 | additive dosage | 0.5 | ... 98 % identity |
| P995 | additive dosage | 0.5 | ... 99.5 % identity |

- Reads: wgsim 2×150 bp, insert 400 bp, 60× over the monoploid genome (15× per homolog), error 0.2 %.
- Replicates: seeds 1–10.

## Methods compared (same reads, same phenotypes)
1. **snp_diploid** – bwa-mem → GATK HaplotypeCaller `-ploidy 2` → GT coded 0/1/2 → GEMMA LMM. (the default many people run)
2. **gwaspoly** – same BAMs → GATK HaplotypeCaller `-ploidy 4` → GT dosage 0–4 → GWASpoly additive model. (the "patch" inside the alignment paradigm)
3. **kmer_pa** – KMERIA k-mer matrix recoded as presence/absence, same LMM (bimbamAsso), same kinship. (what presence/absence k-mer GWAS such as kmersGWAS encodes; this isolates the encoding, it is not the kmersGWAS program)
4. **kmeria** – KMERIA v2.0.4 official steps: count (k=31, canonical) → kctm → filter (depth/ploidy correction) → m2b → bimbamAsso.

## Scoring
- Significance: SNP methods Bonferroni 0.05/M; k-mer methods 0.05/(M/31) (the KMERIA paper's Meff).
- **Detected**: ≥1 significant marker within ±50 kb of the causal site.
- **Ghost**: ≥1 significant marker within ±50 kb of the paralogous copy (signal assigned to the wrong copy).
- **False windows**: significant markers on chr1 outside both windows, counted in 100 kb bins.
- k-mers are placed by exact match against the reference + all founder haplotypes (a perfect pangenome);
  we also report how many significant true-window k-mers do **not** occur in the linear reference.

## Predictions (analytic, n = 400, from the correlation between each encoding and true dosage)
Correlation with true dosage at founder freq 0.5: diploid GT 0.71, presence/absence 0.52, dosage 1.0.

| encoding | expected power at U |
|---|---|
| dosage (gwaspoly, kmeria) | ≈ 95–98 % |
| diploid GT (snp_diploid) | ≈ 35–50 % |
| presence/absence (kmer_pa) | ≈ 5–10 % |

1. **Act 2 (the contradiction)**: at U, snp_diploid and kmer_pa miss a 10 %-variance QTL in a large share of replicates although dosage-aware methods find it nearly always.
2. **Act 3 (patch, then second wall)**: gwaspoly ≈ kmeria at U and P95; gwaspoly loses power at P98 and more at P995 (reads mismap or get MAPQ≈0 between near-identical copies) and may produce ghost hits; kmeria stays ≈ flat across identity because the alt-allele k-mers remain unique.
3. **Act 4 (boundary)**: at CTRL (dominant, freq 0.1) every method detects the QTL; kmer_pa / snp_diploid are expected to be at least as strong as the dosage methods there.
4. Runtime: GATK ploidy 4 genotyping is the most expensive step; kmeria is cheaper end-to-end.

Any outcome that differs from 1–4 is reported as observed.

---

## Amendment 1 (2026-11-04, after the n = 20 smoke test, before any full replicate)
**What the smoke test showed** (n = 20, seed 101, `eval/check_snp_calls.py`): GATK `-ploidy 4` genotypes at the
causal sites correlated with true dosage at r = 1.000 (U, P95, P98) and r = 0.975 (P995); only the site QUAL
dropped (~5.5k vs ~15k). With 2×150 bp reads and 400 bp inserts a fragment spans ~700 bp, so 0.5 % divergence
still leaves ~3–4 distinguishing sites per fragment. Prediction 2 as written is therefore likely to fail at the
planned identity levels.

**Change**: two traits added, everything else unchanged (original predictions 1–4 stand and will be scored):

| trait | QTL model | founder alt freq | where |
|---|---|---|---|
| P999 | additive dosage | 0.5 | paralogous copy at 99.9 % identity (20 differences per 20 kb) |
| P100 | additive dosage | 0.5 | paralogous copy 100 % identical (a recent segmental duplication) |

**Additional predictions**: gwaspoly loses power at P999 and largely fails at P100 (reads from two identical copies
get MAPQ ≈ 0 and are ignored by HaplotypeCaller); kmeria detects both (the causal alt allele exists on one copy only,
so its k-mers stay unique). Replicates: seeds 1–5 first; 6–10 only if needed.
