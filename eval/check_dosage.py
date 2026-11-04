#!/usr/bin/env python3
"""Sanity check: KMERIA copy number of the causal alt-allele k-mers vs the simulated true dosage.
A correlation near 1 also proves that matrix columns follow samples.txt order.
Usage: check_dosage.py <workdir>"""
import glob, sys
import numpy as np

W = sys.argv[1]; K = 31
comp = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(comp)[::-1]
g = np.load(f"{W}/sim/genome.npz")
sites = [l.rstrip("\n").split("\t") for l in open(f"{W}/sim/sites.tsv")][1:]
causal = {s[5]: (s[0], int(s[1]), s[3]) for s in sites if s[4] == "causal"}
true = np.genfromtxt(f"{W}/sim/causal_dosage.tsv", names=True, dtype=None, encoding=None)
spos, chrom, alt, fa = g["spos"], g["chrom"], g["alt"], g["fa"]
want = {}
for name, (c, p, a) in causal.items():
    # alt k-mers on the haplotype of a founder that carries the causal allele
    row = np.flatnonzero((chrom == c) & (spos == p))[0]
    f = np.flatnonzero(fa[row])[0]
    seq = g[f"ref_{c}"].copy(); m = (chrom == c) & fa[:, f]; seq[spos[m]] = alt[m]
    for st in range(p - K + 1, p + 1):
        k = seq[st:st + K].tobytes().decode()
        want[k] = name; want[rc(k)] = name
got = {n: [] for n in causal}
for f in glob.glob(f"{W}/kmer/flt/*.txt"):
    for line in open(f):
        k = line.split("\t", 1)[0][::-1]  # text output is reversed in v2.0.4, see flt_to_pa.py
        if k in want:
            got[want[k]].append(np.array(line.split("\t")[1:], float))
for n in causal:
    if not got[n]:
        print(n, "no alt k-mer survived filtering"); continue
    cn = np.mean(got[n], 0); r = np.corrcoef(cn, true[n])[0, 1]
    print(f"{n:5s} alt k-mers kept={len(got[n]):2d}  r(copy number, true dosage)={r:.3f}")
