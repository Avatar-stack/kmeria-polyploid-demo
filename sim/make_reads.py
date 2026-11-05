#!/usr/bin/env python3
"""Build the P homolog sequences of one individual and simulate paired-end reads with wgsim."""
import argparse, os, pickle, subprocess, tempfile
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("simdir"); ap.add_argument("index", type=int); ap.add_argument("outdir")
ap.add_argument("--depth", type=float, default=60, help="total depth over the monoploid genome")
ap.add_argument("--read-len", type=int, default=150)
ap.add_argument("--err", type=float, default=0.002)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

g = np.load(f"{a.simdir}/genome.npz")
mosaic = pickle.load(open(f"{a.simdir}/mosaic.pkl", "rb"))
chroms = sorted({k[2] for k in mosaic}); P = 1 + max(k[1] for k in mosaic)
spos, chrom, alt, fa = g["spos"], g["chrom"], g["alt"], g["fa"]
s = f"S{a.index:03d}"
tmp = tempfile.mkdtemp(dir=a.outdir)
fa_path, mono = f"{tmp}/{s}.fa", 0
with open(fa_path, "w") as fh:
    for c in chroms:
        m = np.flatnonzero(chrom == c); mono += len(g[f"ref_{c}"])
        for h in range(P):
            starts, fids = mosaic[(a.index, h, c)]
            carry = fa[m, fids[np.searchsorted(starts, spos[m], side="right") - 1]]
            seq = g[f"ref_{c}"].copy(); seq[spos[m][carry]] = alt[m][carry]
            fh.write(f">{c}_h{h}\n{seq.tobytes().decode()}\n")
pairs = int(a.depth * mono / (2 * a.read_len))
r1, r2 = f"{a.outdir}/{s}_1.fq", f"{a.outdir}/{s}_2.fq"
subprocess.run(["wgsim", "-e", str(a.err), "-d", "400", "-s", "40", "-N", str(pairs),
                "-1", str(a.read_len), "-2", str(a.read_len), "-r", "0", "-R", "0", "-X", "0",
                "-S", str(a.seed * 100000 + a.index), fa_path, r1, r2],
               check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.run(["gzip", "-1", "-f", r1, r2], check=True)
os.remove(fa_path); os.rmdir(tmp)
