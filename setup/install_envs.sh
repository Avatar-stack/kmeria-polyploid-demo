#!/usr/bin/env bash
# Create the two environments used by pipelines/run_rep.sh (no sudo needed).
#   kmeriaenv : KMERIA's own kmeria_env.yaml (C++ runtime for the pre-built binaries)
#   gw        : simulation + SNP pipelines (bwa, GATK4, GEMMA, wgsim, R + GWASpoly)
# MIRROR=https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud is much faster from mainland China.
set -euxo pipefail
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}
M=${MICROMAMBA:-$HOME/bin/micromamba}
TOOLS=${TOOLS:-$HOME/kmeria-demo/tools}
MIRROR=${MIRROR:-https://conda.anaconda.org}
if [ ! -x "$M" ]; then mkdir -p "$(dirname "$M")"; curl -sL https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME" bin/micromamba; fi
[ -d "$MAMBA_ROOT_PREFIX/envs/kmeriaenv" ] || "$M" create -y -n kmeriaenv -f "$TOOLS/KMERIA/kmeria_env.yaml"
if [ ! -x "$MAMBA_ROOT_PREFIX/envs/gw/bin/gatk" ]; then
  rm -rf "$MAMBA_ROOT_PREFIX/envs/gw"
  "$M" create -y -n gw --override-channels -c "$MIRROR/conda-forge" -c "$MIRROR/bioconda" \
    python=3.11 numpy pandas scipy matplotlib bwa samtools bcftools gatk4 gemma wgsim htslib \
    r-base=4.4 r-ggplot2 r-matrix r-rrblup r-tidyr r-rlang r-scam
fi
# GWASpoly from a local clone (github.com/jendelman/GWASpoly)
"$M" run -n gw Rscript -e 'if (!requireNamespace("GWASpoly", quietly=TRUE)) install.packages(Sys.getenv("GWASPOLY_SRC"), repos=NULL, type="source")' \
  || true
"$M" run -n gw Rscript -e 'library(GWASpoly); cat("GWASpoly", as.character(packageVersion("GWASpoly")), "\n")'
chmod 755 "$TOOLS"/KMERIA/bin/* "$TOOLS"/KMERIA/external_tools/* "$TOOLS"/KMERIA/bimbamAsso/*
echo ENVS_DONE
