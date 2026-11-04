#!/usr/bin/env python3
"""VCF GT -> alt-allele dosage (0..ploidy). Writes BIMBAM mean-genotype (for GEMMA) and GWASpoly numeric CSV.
Filters: MAF >= 0.05, missing <= 20 %; missing calls are mean-imputed for BIMBAM only."""
import gzip, sys
import numpy as np

vcf, ploidy, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
samples, ids, rows = None, [], []
with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        f = line.rstrip("\n").split("\t")
        if line.startswith("#"):
            samples = f[9:]; continue
        gi = f[8].split(":").index("GT")
        d = []
        for x in f[9:]:
            gt = x.split(":")[gi].replace("|", "/").split("/")
            d.append(np.nan if "." in gt else sum(a != "0" for a in gt))
        d = np.array(d, float)
        miss = np.isnan(d).mean()
        maf = np.nanmean(d) / ploidy if miss < 1 else 0
        if miss <= 0.2 and min(maf, 1 - maf) >= 0.05:
            ids.append((f"{f[0]}:{f[1]}", f[0], f[1], f[3], f[4])); rows.append(d)
G = np.array(rows)
with open(out + ".bimbam", "w") as fh:  # GEMMA mean genotype: scale dosage to 0..2
    for (i, c, p, r, a), d in zip(ids, G):
        d = np.where(np.isnan(d), np.nanmean(d), d) * 2 / ploidy
        fh.write(f"{i},{a},{r}," + ",".join(f"{x:.4g}" for x in d) + "\n")
with open(out + ".gwaspoly.csv", "w") as fh:
    fh.write("Marker,Chrom,Position," + ",".join(samples) + "\n")
    for (i, c, p, r, a), d in zip(ids, G):
        fh.write(f"{i},{c},{p}," + ",".join("NA" if np.isnan(x) else str(int(x)) for x in d) + "\n")
open(out + ".samples", "w").write("\n".join(samples) + "\n")
print(f"{vcf}: kept {len(ids)} markers x {len(samples)} samples", file=sys.stderr)
