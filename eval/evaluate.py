#!/usr/bin/env python3
"""Score one replicate: detection / ghost (signal on the paralog copy) / false windows per trait x method.

k-mer placement: every 31-mer of the reference and of all founder haplotypes is indexed
(founders are colinear with the reference, so this is a perfect "pangenome"). A tested k-mer
gets the set of positions where it (or its reverse complement) occurs.
Thresholds: SNP methods Bonferroni 0.05/M; k-mer methods 0.05/(M/31) as in the KMERIA paper.
"""
import glob, json, os, sys
import numpy as np

W = sys.argv[1]; K = 31; WIN = 50_000; ALPHA = 0.05
truth = json.load(open(f"{W}/sim/truth.json")); seed = truth["args"]["seed"]
CH = {"chr1": 0, "chr2": 2_000_000}  # global offset
LUT = np.full(256, 255, np.uint8)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i


def read_fasta(path):
    name, out = None, {}
    for line in open(path, "rb"):
        if line.startswith(b">"):
            name = line[1:].strip().decode(); out[name] = []
        else:
            out[name].append(line.strip())
    return {k: b"".join(v) for k, v in out.items()}


def codes(seq):  # forward 2-bit codes of all k-mers of a byte string
    a = LUT[np.frombuffer(seq, np.uint8)].astype(np.uint64)
    n = len(a) - K + 1; c = np.zeros(n, np.uint64)
    for j in range(K):
        c = (c << np.uint64(2)) | a[j:j + n]
    return c


def str_codes(kmers):
    b = np.frombuffer("".join(kmers).encode(), np.uint8).reshape(len(kmers), K)
    a = LUT[b].astype(np.uint64)
    fwd = np.zeros(len(kmers), np.uint64); rc = np.zeros(len(kmers), np.uint64)
    for j in range(K):
        fwd = (fwd << np.uint64(2)) | a[:, j]
        rc = (rc << np.uint64(2)) | (np.uint64(3) - a[:, K - 1 - j])
    return fwd, rc


def build_index():
    ref = read_fasta(f"{W}/sim/ref.fa"); fo = read_fasta(f"{W}/sim/founders.fa")
    ref_c = []; all_c = []; all_p = []
    for c, off in CH.items():
        rc = codes(ref[c]); pos = np.arange(len(rc), dtype=np.uint32) + off
        ref_c.append(rc); all_c.append(rc); all_p.append(pos)
        for name, s in fo.items():
            if name.endswith(c):
                fc = codes(s); d = fc != rc  # only k-mers that differ from the reference add information
                all_c.append(fc[d]); all_p.append(pos[d])
    c = np.concatenate(all_c); p = np.concatenate(all_p)
    o = np.lexsort((p, c)); c, p = c[o], p[o]
    keep = np.ones(len(c), bool); keep[1:] = (c[1:] != c[:-1]) | (p[1:] != p[:-1])
    refset = np.unique(np.concatenate(ref_c))
    return c[keep], p[keep], refset


def place(kmers, idx):
    c, p, refset = idx
    fwd, rc = str_codes(kmers)
    first = np.full(len(kmers), -1, np.int64); npos = np.zeros(len(kmers), np.int32)
    for q in (fwd, rc):
        lo = np.searchsorted(c, q, "left"); hi = np.searchsorted(c, q, "right")
        has = hi > lo
        first[has & (first < 0)] = p[lo[has & (first < 0)]]
        npos += (hi - lo).astype(np.int32)
    in_ref = np.isin(fwd, refset) | np.isin(rc, refset)
    return first, npos, in_ref


def load_p(path):
    ids, ps = [], []
    for line in open(path):
        f = line.split()
        if not f or f[0] in ("MARKER_ID", "rs", "marker"):
            continue
        try:
            ps.append(float(f[-1])); ids.append(f[0])
        except ValueError:
            continue
    return ids, np.array(ps)


idx = None
rows, man = [], {}
methods = {"snp_diploid": ("gwas/snp_diploid.{t}.p", "snp"), "gwaspoly": ("gwas/gwaspoly.{t}.p", "snp"),
           "kmer_pa": ("kmer/asso/pa.{t}.ps", "kmer"), "kmeria": ("kmer/asso/dosage.{t}.ps", "kmer")}
for q in truth["qtls"]:
    t = q["trait"]; tw = CH[q["chrom"]] + q["pos"]
    gw = None if q["paralog_pos"] is None else CH[q["chrom"]] + q["paralog_pos"]
    others = [CH[o["chrom"]] + o["pos"] for o in truth["qtls"] if o["trait"] != t]
    for m, (pat, kind) in methods.items():
        f = f"{W}/{pat.format(t=t)}"
        if not os.path.exists(f):
            continue
        ids, ps = load_p(f)
        M = len(ps) if kind == "snp" else len(ps) / K
        thr = ALPHA / M
        if kind == "snp":
            pos = np.array([CH[i.split(":")[0]] + int(i.split(":")[1]) for i in ids]); npos = np.ones(len(ids), int)
            in_ref = np.ones(len(ids), bool)
        else:
            if idx is None:
                idx = build_index()
            pos, npos, in_ref = place(ids, idx)
        sig = ps < thr; placed = pos >= 0
        near = lambda x: placed & (np.abs(pos - x) <= WIN)
        tw_m, gw_m = near(tw), (near(gw) if gw is not None else np.zeros(len(ps), bool))
        on1 = placed & (pos < CH["chr2"])
        fp = sig & on1 & ~tw_m & ~gw_m
        fp_bins = len(np.unique(pos[fp] // (2 * WIN)))
        rows.append(dict(seed=seed, trait=t, identity=q["identity"], model=q["model"], method=m,
                         n_tested=len(ps), threshold=thr,
                         detected=bool((sig & tw_m).any()),
                         min_p_true=float(ps[tw_m].min()) if tw_m.any() else 1.0,
                         n_markers_true=int(tw_m.sum()),
                         ghost=bool((sig & gw_m).any()),
                         min_p_ghost=float(ps[gw_m].min()) if gw_m.any() else None,
                         false_windows=int(fp_bins),
                         false_at_other_qtl=int(sum((sig & near(o)).any() for o in others)),
                         sig_total=int(sig.sum()), sig_unplaced=int((sig & ~placed).sum()),
                         sig_true_multi=int((sig & tw_m & (npos > 1)).sum()),
                         sig_true_not_in_ref=int((sig & tw_m & ~in_ref).sum()),
                         sig_true=int((sig & tw_m).sum())))
        # compact Manhattan data for chr1: keep all markers with p < 1e-2 plus a 2 % sample of the rest
        keep = on1 & ((ps < 1e-2) | (np.random.default_rng(0).random(len(ps)) < 0.02))
        man[f"{t}|{m}"] = dict(pos=pos[keep].tolist(), p=ps[keep].tolist(), thr=thr,
                               multi=(npos[keep] > 1).tolist())
os.makedirs(f"{W}/eval", exist_ok=True)
json.dump(rows, open(f"{W}/eval/metrics.json", "w"), indent=1)
json.dump(man, open(f"{W}/eval/manhattan.json", "w"))
for r in rows:
    print(r["trait"], r["method"], "detected" if r["detected"] else "MISSED",
          f"min_p={r['min_p_true']:.2e} thr={r['threshold']:.1e} ghost={r['ghost']} fp={r['false_windows']}")
