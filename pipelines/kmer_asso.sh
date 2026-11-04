#!/usr/bin/env bash
# k-mer association: kinship from a 0.2 % k-mer sample (bimbamKin), then bimbamAsso per trait for
#   dosage  = KMERIA m2b output (official path)
#   pa      = same filtered k-mers recoded as presence/absence (encoding control)
# Usage: kmer_asso.sh <workdir> "<kmeria env runner>" "<gw env runner>" <repo>
set -euo pipefail
W=$1; KM=$2; GW=$3; REPO=$4; J=${J:-32}
A=$W/kmer/asso; mkdir -p "$A" "$W/kmer/pa"
mapfile -t S < "$W/sim/samples.txt"
printf '%s %s\n' $(for s in "${S[@]}"; do echo "$s $s"; done) > "$A/samples.fam"

# presence/absence bimbam from the filtered (copy-number corrected) text matrix
$GW python "$REPO/pipelines/flt_to_pa.py" "$W/kmer/flt" "$W/kmer/pa"

# kinship from a random 0.2 % of dosage k-mers (KMERIA wiki uses ~0.1 %)
zcat "$W"/kmer/bimbam/*.bimbam.gz | awk 'BEGIN{srand(7)} rand()<0.002' > "$A/sketch.bimbam"
$KM bimbamKin "$A/sketch.bimbam" "$A/kin.txt" -b > "$A/kin.log" 2>&1

NT=$(head -1 "$W/sim/pheno.tsv" | awk -F'\t' '{print NF-1}')
jobs=()
for i in $(seq 1 "$NT"); do
  t=$(head -1 "$W/sim/pheno.tsv" | cut -f$((i+1)))
  awk -F'\t' -v c=$((i+1)) 'NR>1{print $1, $1, $c}' "$W/sim/pheno.tsv" > "$A/$t.pheno"
  for enc in dosage pa; do
    d=$W/kmer/bimbam; [ $enc = pa ] && d=$W/kmer/pa
    for f in "$d"/*.bimbam.gz; do
      jobs+=("$KM bimbamAsso -b -g -t $f -s $A/samples.fam -p $A/$t.pheno -k $A/kin.txt -o $A/tmp.$enc.$t.$(basename "$f" .bimbam.gz)")
    done
  done
done
printf '%s\n' "${jobs[@]}" | xargs -P "$J" -I{} bash -c '{} > /dev/null 2>&1'
for i in $(seq 1 "$NT"); do
  t=$(head -1 "$W/sim/pheno.tsv" | cut -f$((i+1)))
  for enc in dosage pa; do cat "$A"/tmp.$enc.$t.*.ps > "$A/$enc.$t.ps"; done
done
rm -f "$A"/tmp.*
