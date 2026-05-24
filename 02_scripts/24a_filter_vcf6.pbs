#!/bin/bash
#PBS -V
#PBS -N filter_diversity
#PBS -l nodes=node19:ppn=56,mem=64gb,walltime=500:00:00
#PBS -q batch
#PBS -o /usr/scratch2/userdata2/mtofflemire/projects/leuconotopicus/reseq_WD/04_vcf/filter_diversity.log
#PBS -j oe

set -euo pipefail

cd "${PBS_O_WORKDIR:-$PWD}"

source ~/miniconda3/etc/profile.d/conda.sh
conda activate genomics

OUTDIR=/usr/scratch2/userdata2/mtofflemire/projects/leuconotopicus/reseq_WD/04_vcf

SNPVCF=$OUTDIR/leuconotopicus.cohort_raw.SNPs.vcf.gz
OUTVCF=$OUTDIR/leuconotopicus.cohort_filteredQC_no.maf.full.SNPs.vcf.gz

[ -s "$SNPVCF.tbi" ] || tabix -p vcf "$SNPVCF"

bcftools +fill-tags "$SNPVCF" --threads 56 -Ou -- -t AF,F_MISSING \
| bcftools view \
  -m2 \
  -M2 \
  -f .,PASS \
  -e 'AF==1 | AF==0 | F_MISSING > 0.75 | TYPE~"indel" | QUAL < 30.0 | (INFO/DP!="." && (INFO/DP < 12 || INFO/DP > 793)) | (MQ!="." && MQ < 40.0) | (MQRankSum!="." && MQRankSum < -12.5) | (FS!="." && FS > 60.0) | (SOR!="." && SOR > 3.0) | (ReadPosRankSum!="." && ReadPosRankSum < -8.0) | (QD!="." && QD < 2.0)' \
  --threads 56 \
  -Oz \
  -o "$OUTVCF"

tabix -f -p vcf "$OUTVCF"

printf "Filtered SNP count: "
bcftools index -n "$OUTVCF"

echo
echo "Finished:"
echo "$OUTVCF"