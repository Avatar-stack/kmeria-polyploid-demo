#!/usr/bin/env bash
# GEMMA LMM (Wald) for every trait column. Usage: gemma_bimbam.sh geno.bimbam pheno.tsv outprefix "<env runner>"
set -euo pipefail
G=$1; PH=$2; OUT=$3; RUN=${4:-}
D=$(dirname "$OUT"); B=$(basename "$OUT"); cd "$D"
NT=$(head -1 "$PH" | awk -F'\t' '{print NF-1}')
# phenotype rows must follow the genotype sample order (written next to the bimbam by vcf_to_dosage.py)
awk -F'\t' 'NR==FNR{if(FNR>1){$1=$1; v[$1]=$0}; next} {split(v[$1],a," "); s=""; for(i=2;i<=length(a);i++) s=s (i>2?" ":"") a[i]; print s}' \
  OFS=' ' "$PH" "${G%.bimbam}.samples" > "$B.pheno"
$RUN gemma -g "$G" -p "$B.pheno" -gk 1 -o "$B.kin" -outdir . > /dev/null
for i in $(seq 1 "$NT"); do
  t=$(head -1 "$PH" | cut -f$((i+1)))
  $RUN gemma -g "$G" -p "$B.pheno" -n "$i" -k "$B.kin.cXX.txt" -lmm 1 -o "$B.$t" -outdir . > /dev/null
  awk 'NR>1{print $2"\t"$NF}' "$B.$t.assoc.txt" > "$B.$t.p"   # marker id, p_wald
done
