#!/usr/bin/env python3
"""Simulate an autotetraploid GWAS population with planted QTLs inside paralogous segments.

Design (see PREREGISTRATION.md):
  chr1 (2 Mb): one dosage QTL in unique sequence, one dominant control QTL, and three dosage
               QTLs inside 20 kb segments that have a paralogous copy elsewhere on chr1 at
               95% / 98% / 99.5% identity (the single "knob").
  chr2 (0.5 Mb): polygenic background only.
Only substitutions are simulated, so every founder haplotype is colinear with the reference:
a k-mer found in a founder has an exact reference coordinate (our stand-in for a pangenome).
"""
import argparse, json, os, pickle
import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
CHROMS = {"chr1": 2_000_000, "chr2": 500_000}
SEG = 20_000
# name, model, chrom, causal pos, paralog copy start (None = unique), identity
QTLS = [
    ("U",    "dosage",   "chr1",   250_000, None,      1.0),
    ("CTRL", "dominant", "chr1",   450_000, None,      1.0),
    ("P95",  "dosage",   "chr1",   710_000, 1_500_000, 0.95),
    ("P98",  "dosage",   "chr1",   910_000, 1_650_000, 0.98),
    ("P995", "dosage",   "chr1", 1_110_000, 1_800_000, 0.995),
]


def random_seq(rng, n, gc=0.36):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return BASES[rng.choice(4, size=n, p=p)]


def other_base(rng, ref, avoid=None):
    choices = [b for b in BASES if b != ref and b != avoid]
    return choices[rng.integers(len(choices))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-ind", type=int, default=400)
    ap.add_argument("--n-founders", type=int, default=60)
    ap.add_argument("--ploidy", type=int, default=4)
    ap.add_argument("--snp-rate", type=float, default=1 / 150)
    ap.add_argument("--seg-mean", type=float, default=25_000, help="mean mosaic block length (bp)")
    ap.add_argument("--h2-major", type=float, default=0.10)
    ap.add_argument("--h2-bg", type=float, default=0.20)
    ap.add_argument("--n-bg", type=int, default=100)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.out, exist_ok=True)
    F, P, N = a.n_founders, a.ploidy, a.n_ind

    # 1. ancestral reference with paralogous copies
    ref = {c: random_seq(rng, L) for c, L in CHROMS.items()}
    para = {}
    for name, model, c, pos, copy, ident in QTLS:
        if copy is None:
            continue
        src = (pos // SEG) * SEG
        seg = ref[c][src:src + SEG].copy()
        idx = rng.choice(SEG, size=int(round(SEG * (1 - ident))), replace=False)
        for i in idx:
            seg[i] = other_base(rng, seg[i])
        ref[c][copy:copy + SEG] = seg
        para[name] = dict(src=src, copy=copy, n_diff=len(idx))

    # 2. polymorphic sites (+ forced causal sites)
    causal = {(c, pos): (name, copy) for name, model, c, pos, copy, _ in QTLS}
    rows = []  # chrom, pos, ref, alt, kind, name
    for c, L in CHROMS.items():
        for p in np.flatnonzero(rng.random(L) < a.snp_rate):
            if (c, int(p)) not in causal:
                rows.append((c, int(p), ref[c][p], other_base(rng, ref[c][p]), "snp", ""))
    for (c, p), (name, copy) in causal.items():
        # alt must also differ from the paralog base, so alt-allele k-mers stay unique
        avoid = None if copy is None else ref[c][copy + (p - (p // SEG) * SEG)]
        rows.append((c, p, ref[c][p], other_base(rng, ref[c][p], avoid), "causal", name))
    rows.sort(key=lambda r: (r[0], r[1]))
    S = len(rows)
    chrom = np.array([r[0] for r in rows]); spos = np.array([r[1] for r in rows])
    alt = np.array([r[3] for r in rows], dtype=np.uint8)
    kind = np.array([r[4] for r in rows])

    target = {"U": 0.5, "P95": 0.5, "P98": 0.5, "P995": 0.5, "CTRL": 0.1}
    fa = np.zeros((S, F), dtype=bool)  # founder carries alt
    for i, r in enumerate(rows):
        f = target[r[5]] if r[4] == "causal" else np.clip(rng.beta(0.4, 0.4), 2 / F, 1 - 2 / F)
        fa[i, rng.choice(F, size=int(round(f * F)), replace=False)] = True

    # 3. individuals: each homolog of each chromosome is a mosaic of founder blocks
    mosaic = {}  # (ind, hom, chrom) -> (block starts, founder ids)
    dos = np.zeros((S, N), dtype=np.int8)
    for c, L in CHROMS.items():
        m = np.flatnonzero(chrom == c)
        for i in range(N):
            for h in range(P):
                starts = [0]
                while (nxt := starts[-1] + int(rng.exponential(a.seg_mean)) + 1) < L:
                    starts.append(nxt)
                starts = np.array(starts); fids = rng.integers(F, size=len(starts))
                mosaic[(i, h, c)] = (starts, fids)
                dos[m, i] += fa[m, fids[np.searchsorted(starts, spos[m], side="right") - 1]]

    # 4. phenotypes: one trait per QTL = major QTL + shared chr2 background + noise
    def scale(x, v):
        x = x - x.mean()
        return x / x.std() * np.sqrt(v)
    freq = fa.mean(1)
    bg = rng.choice(np.flatnonzero((chrom == "chr2") & (kind == "snp") & (freq > 0.1) & (freq < 0.9)),
                    size=a.n_bg, replace=False)
    g_bg = scale(rng.normal(size=a.n_bg) @ dos[bg].astype(float), a.h2_bg)
    name2row = {r[5]: i for i, r in enumerate(rows) if r[4] == "causal"}
    pheno, truth = {}, []
    for name, model, c, pos, copy, ident in QTLS:
        d = dos[name2row[name]].astype(float)
        g = d if model == "dosage" else (d > 0).astype(float)
        pheno[name] = scale(g, a.h2_major) + g_bg + rng.normal(scale=np.sqrt(1 - a.h2_major - a.h2_bg), size=N)
        truth.append(dict(trait=name, model=model, chrom=c, pos=pos, identity=ident,
                          paralog_pos=None if copy is None else copy + (pos - (pos // SEG) * SEG),
                          alt_freq_founders=float(freq[name2row[name]]),
                          dose_counts=np.bincount(dos[name2row[name]], minlength=P + 1).tolist()))
    # self-check: realised causal dosage must not be degenerate
    for t in truth:
        assert max(t["dose_counts"]) < N, t

    # 5. outputs
    samples = [f"S{i:03d}" for i in range(N)]
    with open(f"{a.out}/ref.fa", "w") as fh:
        for c in CHROMS:
            fh.write(f">{c}\n{ref[c].tobytes().decode()}\n")
    with open(f"{a.out}/sites.tsv", "w") as fh:
        fh.write("chrom\tpos\tref\talt\tkind\tname\tfounder_freq\n")
        for i, r in enumerate(rows):
            fh.write(f"{r[0]}\t{r[1]}\t{chr(r[2])}\t{chr(r[3])}\t{r[4]}\t{r[5]}\t{freq[i]:.4f}\n")
    with open(f"{a.out}/pheno.tsv", "w") as fh:
        fh.write("sample\t" + "\t".join(pheno) + "\n")
        for i, s in enumerate(samples):
            fh.write(s + "\t" + "\t".join(f"{pheno[t][i]:.6f}" for t in pheno) + "\n")
    with open(f"{a.out}/causal_dosage.tsv", "w") as fh:
        fh.write("sample\t" + "\t".join(t["trait"] for t in truth) + "\n")
        for i, s in enumerate(samples):
            fh.write(s + "\t" + "\t".join(str(dos[name2row[t['trait']], i]) for t in truth) + "\n")
    open(f"{a.out}/samples.txt", "w").write("\n".join(samples) + "\n")
    json.dump(dict(args=vars(a), qtls=truth, n_sites=S, paralogs=para, bg_rows=bg.tolist()),
              open(f"{a.out}/truth.json", "w"), indent=1)
    np.savez_compressed(f"{a.out}/genome.npz", spos=spos, chrom=chrom, alt=alt, fa=fa,
                        **{f"ref_{c}": ref[c] for c in CHROMS})
    pickle.dump(mosaic, open(f"{a.out}/mosaic.pkl", "wb"))
    with open(f"{a.out}/founders.fa", "w") as fh:  # our "pangenome" for placing k-mers
        for f in range(F):
            for c in CHROMS:
                seq = ref[c].copy(); m = (chrom == c) & fa[:, f]
                seq[spos[m]] = alt[m]
                fh.write(f">F{f:02d}_{c}\n{seq.tobytes().decode()}\n")
    print(json.dumps([{k: t[k] for k in ("trait", "dose_counts")} for t in truth]))


if __name__ == "__main__":
    main()
