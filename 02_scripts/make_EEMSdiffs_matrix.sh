#need to make genomic dissimilarity matrix with the 1-ibs command in PLINK1.9 
plink \
  --vcf '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data/Leuconotopicus_albolarvatus_pruned_snps.vcf.gz' \
  --double-id \
  --allow-extra-chr \
  --distance square 1-ibs \
  --out '/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/02_scripts/EEMS/eems.Leuconotopicus.diffs'


