#!/bin/bash

############################################################
## Paths
############################################################

META="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data/Leuconotopicus_albolarvatus_Meta.csv"

VCF="/usr/scratch2/userdata2/mtofflemire/projects/leuconotopicus/reseq_WD/04_vcf/leuconotopicus.cohort_filteredQC_no.maf.full.SNPs.vcf.gz"

BASE="/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs"

POPDIR="$BASE/02_scripts/pop_files"

OUTDIR="$BASE/03_results/05_genome_wide_variation"

PIDIR="$OUTDIR/pi"
TAJDIR="$OUTDIR/tajimasD"
FSTDIR="$OUTDIR/pairwise_fst"

mkdir -p "$POPDIR"
mkdir -p "$PIDIR"
mkdir -p "$TAJDIR"
mkdir -p "$FSTDIR"

############################################################
## Create population files
############################################################

awk -F',' '

NR==1 {

for(i=1;i<=NF;i++) {
if($i=="SequenceID") seq=i
if($i=="Ecoregion") eco=i

}

next

}

$seq!="" && $eco!="" {

eco_name=$eco

gsub(/[^A-Za-z0-9_]/,"_",eco_name)

print $seq >> "'"$POPDIR"'/" eco_name ".txt"

}

' "$META"

############################################################
## PI + Tajima D
############################################################

for POP in "$POPDIR"/*.txt
do

NAME=$(basename "$POP" .txt)

echo
echo "Running PI: $NAME"

vcftools \
--gzvcf "$VCF" \
--keep "$POP" \
--window-pi 100000 \
--window-pi-step 10000 \
--out "$PIDIR/$NAME"

echo
echo "Running TajimaD: $NAME"

vcftools \
--gzvcf "$VCF" \
--keep "$POP" \
--TajimaD 100000 \
--out "$TAJDIR/$NAME"

done

############################################################
## Pairwise FST
############################################################

POPS=("$POPDIR"/*.txt)

for ((i=0;i<${#POPS[@]};i++))
do

for ((j=i+1;j<${#POPS[@]};j++))
do

POP1=${POPS[i]}
POP2=${POPS[j]}

NAME1=$(basename "$POP1" .txt)
NAME2=$(basename "$POP2" .txt)

echo
echo "Running FST: $NAME1 vs $NAME2"

vcftools \
--gzvcf "$VCF" \
--weir-fst-pop "$POP1" \
--weir-fst-pop "$POP2" \
--fst-window-size 100000 \
--fst-window-step 10000 \
--out "$FSTDIR/${NAME1}_vs_${NAME2}"

done

done

echo
echo "Finished all analyses"