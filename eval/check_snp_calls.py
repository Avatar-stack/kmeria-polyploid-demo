#!/usr/bin/env python3
"""Sanity check: GATK genotype at each causal site vs true dosage (ploidy 2 and 4 calls).
Usage: check_snp_calls.py <workdir>"""
import gzip, json, sys
import numpy as np

W = sys.argv[1]
truth = json.load(open(f"{W}/sim/truth.json"))
true = np.genfromtxt(f"{W}/sim/causal_dosage.tsv", names=True, dtype=None, encoding=None)
order = {s: i for i, s in enumerate(true["sample"])}
for pl in (2, 4):
    want = {(q["chrom"], q["pos"] + 1): q["trait"] for q in truth["qtls"]}  # VCF is 1-based
    found = {}
    with gzip.open(f"{W}/vcf/raw{pl}.vcf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]; continue
            key = (f[0], int(f[1]))
            if key in want:
                gi = f[8].split(":").index("GT")
                d = [np.nan if "." in x.split(":")[gi] else sum(a != "0" for a in x.split(":")[gi].replace("|", "/").split("/")) for x in f[9:]]
                found[want[key]] = (np.array(d, float), f[5])
    for q in truth["qtls"]:
        t = q["trait"]
        if t not in found:
            print(f"ploidy={pl} {t:5s} causal SNP NOT CALLED"); continue
        d, qual = found[t]; td = np.array([true[t][order[s]] for s in samples], float)
        ok = ~np.isnan(d)
        r = np.corrcoef(d[ok], td[ok])[0, 1] if ok.sum() > 2 and d[ok].std() > 0 else float("nan")
        print(f"ploidy={pl} {t:5s} QUAL={qual:>8s} called={ok.mean():.0%} r(GT dosage, true)={r:.3f}")
