#!/bin/bash

set -euo pipefail

############################################################
## Paths
############################################################

VCF="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data/leuconotopicus.cohort_filteredQC_no.maf.full.SNPs.vcf.gz"

BASE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs"

POPDIR="$BASE/02_scripts"

OUTDIR="$BASE/03_results/05_genome_wide_variation_subspecies"

PIDIR="$OUTDIR/pi"

TAJDIR="$OUTDIR/tajimasD"

FSTDIR="$OUTDIR/pairwise_fst"

mkdir -p "$PIDIR"
mkdir -p "$TAJDIR"
mkdir -p "$FSTDIR"

############################################################
## Population files
############################################################

POPS=(

"$POPDIR/albolarvatus_CA_samps.txt"

"$POPDIR/albolarvatus_WA_OR_samps.txt"

"$POPDIR/gravirostris_samps.txt"

)

############################################################
## PI + TajimaD
############################################################

for POP in "${POPS[@]}"
do

NAME=$(basename "$POP" .txt)

echo
echo "Running PI: $NAME"

vcftools \
--gzvcf "$VCF" \
--keep "$POP" \
--window-pi 50000 \
--window-pi-step 50000 \
--out "$PIDIR/$NAME"

echo
echo "Running TajimaD: $NAME"

vcftools \
--gzvcf "$VCF" \
--keep "$POP" \
--TajimaD 50000 \
--out "$TAJDIR/$NAME"

done

############################################################
## Pairwise FST
############################################################

for ((i=0;i<${#POPS[@]};i++))
do

for ((j=i+1;j<${#POPS[@]};j++))
do

POP1="${POPS[i]}"
POP2="${POPS[j]}"

NAME1=$(basename "$POP1" .txt)

NAME2=$(basename "$POP2" .txt)

echo
echo "Running FST: $NAME1 vs $NAME2"

vcftools \
--gzvcf "$VCF" \
--weir-fst-pop "$POP1" \
--weir-fst-pop "$POP2" \
--fst-window-size 50000 \
--fst-window-step 50000 \
--out "$FSTDIR/${NAME1}_vs_${NAME2}"

done

done

echo
echo "Finished all subspecies analyses"