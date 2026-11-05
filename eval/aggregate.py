#!/usr/bin/env python3
"""Aggregate replicates -> results/*.tsv and figures/*.png.
Usage: aggregate.py <runs_dir> <repo_dir> [manhattan_seed]"""
import glob, json, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R, REPO = sys.argv[1], sys.argv[2]
MSEED = sys.argv[3] if len(sys.argv) > 3 else None
os.makedirs(f"{REPO}/results", exist_ok=True); os.makedirs(f"{REPO}/figures", exist_ok=True)
METH = ["snp_diploid", "gwaspoly", "kmer_pa", "kmeria"]
LAB = {"snp_diploid": "SNP, diploid calls (GATK ploidy 2 + GEMMA)", "gwaspoly": "SNP, tetraploid dosage (GATK ploidy 4 + GWASpoly)",
       "kmer_pa": "k-mer presence/absence (same matrix, same LMM)", "kmeria": "KMERIA (k-mer copy number)"}
COL = {"snp_diploid": "#9aa0a6", "gwaspoly": "#1a73e8", "kmer_pa": "#f29900", "kmeria": "#d93025"}
ORDER = ["U", "P95", "P98", "P995", "P999", "P100"]
XL = {"U": "unique", "P95": "95%", "P98": "98%", "P995": "99.5%", "P999": "99.9%", "P100": "100%"}

seeds = sorted(d for d in glob.glob(f"{R}/seed*") if os.path.exists(f"{d}/eval/metrics.json"))
m = pd.DataFrame([r for d in seeds for r in json.load(open(f"{d}/eval/metrics.json"))])
m["neglog10p"] = -np.log10(m["min_p_true"].clip(lower=1e-300))
m["neglog10thr"] = -np.log10(m["threshold"])
summ = (m.groupby(["trait", "method"])
          .agg(reps=("seed", "nunique"), power=("detected", "mean"), ghost_rate=("ghost", "mean"),
               false_windows=("false_windows", "mean"), median_neglog10p=("neglog10p", "median"),
               neglog10thr=("neglog10thr", "mean"), n_tested=("n_tested", "mean"),
               sig_true_not_in_ref=("sig_true_not_in_ref", "sum"), sig_true=("sig_true", "sum"))
          .reset_index())
summ.to_csv(f"{REPO}/results/summary.tsv", sep="\t", index=False, float_format="%.4g")
m.to_csv(f"{REPO}/results/metrics_all.tsv", sep="\t", index=False, float_format="%.4g")

# runtime per stage (wall seconds), mean over replicates
t = []
for d in seeds:
    for line in open(f"{d}/logs/time.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) == 3 and not f[0].startswith("Command"):
            t.append(dict(seed=d, stage=f[0], sec=float(f[1]), maxrss_gb=float(f[2]) / 1e6))
t = pd.DataFrame(t).groupby(["seed", "stage"]).sum().reset_index().groupby("stage").mean(numeric_only=True)
GROUP = {"bwa": "SNP", "hc2": "SNP", "gt2": "SNP", "hc4": "SNP", "gt4": "SNP", "dosage": "SNP", "gemma": "SNP",
         "gwaspoly": "SNP", "kcount": "k-mer", "kctm": "k-mer", "kfilter": "k-mer", "m2b": "k-mer", "kasso": "k-mer"}
t["route"] = [GROUP.get(s, "shared") for s in t.index]
t.to_csv(f"{REPO}/results/runtime.tsv", sep="\t", float_format="%.1f")

# Figure 1: power vs paralog identity
fig, ax = plt.subplots(figsize=(7.5, 4.2))
x = np.arange(len(ORDER))
for me in METH:
    s = summ[summ.method == me].set_index("trait").reindex(ORDER)
    ax.plot(x, s.power, "-o", color=COL[me], label=LAB[me], lw=2.2 if me == "kmeria" else 1.6)
ax.set_xticks(x, [XL[o] for o in ORDER]); ax.set_ylim(-0.04, 1.04)
ax.set_xlabel("where the dosage QTL sits (identity to its paralogous copy)"); ax.set_ylabel(f"detection rate ({len(seeds)} replicates)")
ax.legend(fontsize=8, loc="lower left", frameon=False); ax.grid(alpha=.25)
fig.tight_layout(); fig.savefig(f"{REPO}/figures/fig1_power_vs_identity.png", dpi=200)

# Figure 2: signal strength at the true locus relative to the threshold (per replicate)
fig, axes = plt.subplots(1, len(ORDER) + 1, figsize=(13, 3.6), sharey=True)
for ax, tr in zip(axes, ORDER + ["CTRL"]):
    for j, me in enumerate(METH):
        s = m[(m.trait == tr) & (m.method == me)]
        ax.scatter(np.full(len(s), j) + np.random.default_rng(j).uniform(-.12, .12, len(s)),
                   s.neglog10p - s.neglog10thr, color=COL[me], s=18)
    ax.axhline(0, color="k", lw=.8, ls="--"); ax.set_title(XL.get(tr, "dominant control"), fontsize=9)
    ax.set_xticks(range(4), ["dip", "poly", "P/A", "KMERIA"], fontsize=7, rotation=45)
axes[0].set_ylabel("-log10(p at true locus) - (-log10 threshold)")
fig.tight_layout(); fig.savefig(f"{REPO}/figures/fig2_signal_margin.png", dpi=200)

# Figure 3: runtime by route
fig, ax = plt.subplots(figsize=(6, 3.2))
tt = t[t.route != "shared"].sort_values("sec")
ax.barh(tt.index, tt.sec / 60, color=["#1a73e8" if r == "SNP" else "#d93025" for r in tt.route])
ax.set_xlabel("wall-clock minutes per replicate (n = 400, 56 parallel jobs)")
tot = t.groupby("route").sec.sum() / 60
ax.set_title(f"SNP route {tot.get('SNP', 0):.0f} min vs k-mer route {tot.get('k-mer', 0):.0f} min", fontsize=10)
fig.tight_layout(); fig.savefig(f"{REPO}/figures/fig3_runtime.png", dpi=200)

# Figure 4: Manhattan (chr1) for one replicate, trait P100
if MSEED:
    man = json.load(open(f"{R}/seed{MSEED}/eval/manhattan.json"))
    tr = json.load(open(f"{R}/seed{MSEED}/sim/truth.json"))
    q = {x["trait"]: x for x in tr["qtls"]}
    for trait in ("U", "P100"):
        fig, axes = plt.subplots(4, 1, figsize=(9, 7), sharex=True)
        for ax, me in zip(axes, METH):
            d = man.get(f"{trait}|{me}")
            if d is None:
                continue
            pos = np.array(d["pos"]) / 1e6; lp = -np.log10(np.clip(d["p"], 1e-300, 1))
            ax.scatter(pos, lp, s=3, color=COL[me], rasterized=True)
            ax.axhline(-np.log10(d["thr"]), color="k", lw=.7, ls="--")
            ax.axvline(q[trait]["pos"] / 1e6, color="green", lw=.8, alpha=.6)
            if q[trait]["paralog_pos"]:
                ax.axvline(q[trait]["paralog_pos"] / 1e6, color="purple", lw=.8, ls=":", alpha=.8)
            ax.set_ylabel("-log10 p", fontsize=8); ax.set_title(LAB[me], fontsize=8, loc="left")
        axes[-1].set_xlabel("chr1 position (Mb); green = causal site, purple dotted = paralogous copy")
        fig.tight_layout(); fig.savefig(f"{REPO}/figures/fig4_manhattan_{trait}_seed{MSEED}.png", dpi=170)
print(summ.pivot(index="trait", columns="method", values="power").reindex(ORDER + ["CTRL"])[METH].to_string())
print(t[["sec", "route"]].to_string())
