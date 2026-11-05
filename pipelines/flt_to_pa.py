#!/usr/bin/env python3
"""Recode KMERIA's filtered copy-number matrix as presence/absence (count > 0), BIMBAM gz.
Keeps only k-mers that are present in >= 5 and absent in >= 5 samples (presence/absence k-mer GWAS
tests only variable k-mers). Usage: flt_to_pa.py <flt_dir> <out_dir>"""
import glob, gzip, os, sys
import numpy as np
import pandas as pd

src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
for f in sorted(glob.glob(f"{src}/*.txt")):
    kept = total = 0
    with gzip.open(f"{out}/{os.path.basename(f)[:-4]}.bimbam.gz", "wt", compresslevel=1) as fh:
        for chunk in pd.read_csv(f, sep="\t", header=None, chunksize=200_000):
            # KMERIA v2.0.4 `filter --output-format text` prints k-mers reversed (binary_to_string
            # decodes low bits first); reverse back so ids match the m2b/compressed path.
            ids = np.array([k[::-1] for k in chunk.iloc[:, 0]]); pa = (chunk.iloc[:, 1:].to_numpy() > 0).astype(np.int8)
            n_p = pa.sum(1); ok = (n_p >= 5) & (pa.shape[1] - n_p >= 5)
            total += len(ids); kept += int(ok.sum())
            for k, row in zip(ids[ok], pa[ok]):
                fh.write(k + ",X,Y," + ",".join("2" if x else "0" for x in row) + "\n")
    print(f"{f}: {kept}/{total} variable presence/absence k-mers", file=sys.stderr)
