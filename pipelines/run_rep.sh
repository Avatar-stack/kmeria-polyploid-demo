#!/usr/bin/env bash
# One replicate: simulate -> reads -> 4 GWAS pipelines -> evaluation.
# Usage: pipelines/run_rep.sh <seed> <workdir> [n_ind] [jobs]
# Each stage writes a .done marker, so a failed run can be restarted.
set -euo pipefail
SEED=$1; W=$(realpath -m "$2"); N=${3:-400}; J=${4:-40}
REPO=$(cd "$(dirname "$0")/.." && pwd)
DEPTH=60; RL=150; K=31; PLOIDY=4
ENVS=${MAMBA_ROOT_PREFIX:-$HOME/micromamba}/envs
KMERIA=$REPO/../tools/KMERIA   # repo is cloned next to tools/ on the VM
[ -d "$KMERIA" ] || KMERIA=$HOME/kmeria-demo/tools/KMERIA
# call env binaries directly: `micromamba run` takes a global lock and serialises parallel jobs
GW="env PATH=$ENVS/gw/bin:$PATH"
KM="env LD_LIBRARY_PATH=$KMERIA/lib:$ENVS/kmeriaenv/lib PATH=$KMERIA/bin:$KMERIA/bimbamAsso:$PATH"
mkdir -p "$W"/{sim,reads,bam,gvcf2,gvcf4,vcf,gwas,kmer,logs}
T() { /usr/bin/time -f "$1\t%e\t%M" -a -o "$W/logs/time.tsv" "${@:2}"; }  # stage wall-s, maxRSS-KB
done_() { [ -f "$W/.$1.done" ]; }; mark() { touch "$W/.$1.done"; }
log() { echo "[$(date +%T)] seed=$SEED $*"; }

if ! done_ sim; then log simulate
  T simulate $GW python "$REPO/sim/simulate.py" --out "$W/sim" --seed "$SEED" --n-ind "$N" > "$W/logs/sim.json"
  mark sim; fi

if ! done_ reads; then log reads
  T reads bash -c "seq 0 $((N-1)) | xargs -P $J -I{} $GW python '$REPO/sim/make_reads.py' '$W/sim' {} '$W/reads' --depth $DEPTH --seed $SEED"
  mark reads; fi
mapfile -t S < "$W/sim/samples.txt"

# ---------- SNP pipelines: bwa + GATK HaplotypeCaller (ploidy 2 and 4) ----------
if ! done_ bam; then log bwa
  $GW bwa index "$W/sim/ref.fa" 2>/dev/null; $GW samtools faidx "$W/sim/ref.fa"
  $GW gatk CreateSequenceDictionary -R "$W/sim/ref.fa" >/dev/null 2>&1
  T bwa bash -c "printf '%s\n' ${S[*]} | xargs -P $((J/4)) -I{} bash -c '$GW bwa mem -t 4 -R \"@RG\tID:{}\tSM:{}\" $W/sim/ref.fa $W/reads/{}_1.fq.gz $W/reads/{}_2.fq.gz 2>/dev/null | $GW samtools sort -@2 -o $W/bam/{}.bam - && $GW samtools index $W/bam/{}.bam'"
  mark bam; fi

for PL in 2 4; do
  if ! done_ hc$PL; then log "HaplotypeCaller ploidy=$PL"
    T hc$PL bash -c "printf '%s\n' ${S[*]} | xargs -P $J -I{} $GW gatk --java-options -Xmx2g HaplotypeCaller -R $W/sim/ref.fa -I $W/bam/{}.bam -ploidy $PL -ERC GVCF -O $W/gvcf$PL/{}.g.vcf.gz --native-pair-hmm-threads 1 >/dev/null 2>&1"
    mark hc$PL; fi
  if ! done_ gt$PL; then log "joint genotyping ploidy=$PL"
    printf '%s\n' "${S[@]}" | awk -v d="$W/gvcf$PL" '{print $1"\t"d"/"$1".g.vcf.gz"}' > "$W/vcf/map$PL.txt"
    rm -rf "$W/vcf/db$PL"
    T gt$PL bash -c "$GW gatk --java-options -Xmx16g GenomicsDBImport --sample-name-map $W/vcf/map$PL.txt --genomicsdb-workspace-path $W/vcf/db$PL -L chr1 -L chr2 --batch-size 50 --reader-threads 4 >/dev/null 2>&1 && $GW gatk --java-options -Xmx16g GenotypeGVCFs -R $W/sim/ref.fa -V gendb://$W/vcf/db$PL -O $W/vcf/raw$PL.vcf.gz >/dev/null 2>&1"
    $GW bcftools view -m2 -M2 -v snps -i 'QUAL>=30' "$W/vcf/raw$PL.vcf.gz" -Oz -o "$W/vcf/snp$PL.vcf.gz"
    mark gt$PL; fi
done

if ! done_ snpgwas; then log "SNP GWAS (diploid GT -> GEMMA; ploidy4 GT -> GWASpoly)"
  T dosage $GW python "$REPO/pipelines/vcf_to_dosage.py" "$W/vcf/snp2.vcf.gz" 2 "$W/gwas/dip"
  T dosage $GW python "$REPO/pipelines/vcf_to_dosage.py" "$W/vcf/snp4.vcf.gz" 4 "$W/gwas/tet"
  T gemma bash "$REPO/pipelines/gemma_bimbam.sh" "$W/gwas/dip.bimbam" "$W/sim/pheno.tsv" "$W/gwas/snp_diploid" "$GW"
  T gwaspoly $GW Rscript "$REPO/pipelines/gwaspoly.R" "$W/gwas/tet.gwaspoly.csv" "$W/sim/pheno.tsv" "$W/gwas/gwaspoly" 8
  mark snpgwas; fi

# ---------- k-mer pipelines: KMERIA (dosage) and the same matrix as presence/absence ----------
if ! done_ kcount; then log "kmeria count"
  T kcount bash -c "printf '%s\n' ${S[*]} | xargs -P $((J/4)) -I{} $KM kmeria count -k $K -t 4 -m 3 -M 1000 -o $W/kmer/{}.bin $W/reads/{}_1.fq.gz $W/reads/{}_2.fq.gz >/dev/null 2>&1"
  mark kcount; fi
if ! done_ kmatrix; then log "kctm + filter + m2b"
  for s in "${S[@]}"; do echo "$W/kmer/$s.bin"; done > "$W/kmer/db.list"
  KD=$(python3 -c "print($DEPTH*($RL-$K+1)/$RL)")    # k-mer depth, README FAQ (2)
  for s in "${S[@]}"; do printf '%s\t%s\t%s\n' "$s" "$KD" "$PLOIDY"; done > "$W/kmer/depth.txt"
  rm -rf "$W/kmer/mat" "$W/kmer/flt" "$W/kmer/bimbam"; mkdir -p "$W/kmer/mat"
  T kctm $KM kmeria kctm -i "$W/kmer/db.list" -o "$W/kmer/mat/matrix" -j 16 > "$W/logs/kctm.log" 2>&1
  # official path (as in kmeria_wrapper.pl): compressed filter output -> m2b
  T kfilter $KM kmeria filter -i "$W/kmer/mat" -o "$W/kmer/flt_c" -d "$W/kmer/depth.txt" -p $PLOIDY -t 16 --output-format compressed > "$W/logs/filter.log" 2>&1
  T m2b $KM kmeria m2b --in "$W/kmer/flt_c" --out "$W/kmer/bimbam" --threads 16 > "$W/logs/m2b.log" 2>&1
  # text copy of the same filtered copy numbers, only used to build the presence/absence control.
  # NOTE (v2.0.4): text output writes k-mer strings reversed; flt_to_pa.py reverses them back.
  T kfilter_txt $KM kmeria filter -i "$W/kmer/mat" -o "$W/kmer/flt" -d "$W/kmer/depth.txt" -p $PLOIDY -t 16 --output-format text >> "$W/logs/filter.log" 2>&1
  mark kmatrix; fi
if ! done_ kasso; then log "kinship + association (dosage and presence/absence)"
  T kasso bash "$REPO/pipelines/kmer_asso.sh" "$W" "$KM" "$GW" "$REPO"
  mark kasso; fi

if ! done_ eval; then log evaluate
  T eval $GW python "$REPO/eval/evaluate.py" "$W" > "$W/logs/eval.log"
  mark eval; fi
log finished
