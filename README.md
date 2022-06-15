# variant_viewer
Code for the development of a pathogenic variant viewer app
Variant viewer app displays table of possibly pathogenic variants and an IGV instance as well as some additional information. This is organized in a R markdown html document with three Rshiny apps.
  1. Displays an IGV instance with a pre loaded bed file showing the postions of variants in the proband and Parents
  2. A table of possibly pathogenic variants and their characteristics 
  3. A data dictionary, summary of VEP annotation, and an overview of the process to create the app
  
Inputs for the app are created from the following processes.
  1. VCF files for the family are merged.
  2. VCF files are filtered using vcf and bcftools. Remove sites with missing genotypes for Proband or low quality. (Also need to remove sites where proband is heterogenoeous and parents have no variant
  3. Annotate VCF using VEP (offline no HGVS)
  4. Filter VEP annotation towards sites more likely to be pathogenic
  5. Get trimmed VCF file from postions vcftools --vcf inputs/allfam.vcf --positions inputs/ft_positions.txt --recode --out outputs/final_filter_allfam
  6. annotate trimmed VCF file with HGVS and other data not availble to offline cache.

Table is obtained through annotation and filtering. Filtering scheme is present in annotations.r

## VEP annotation

annotation was command line


ensembl-vep/./vep -i ensembl-vep/inputs/allfam.vcf -o vep_anno.tsv --cache  --everything --dir_cache $HOME/cache "exact" --tab

ensembl-vep/./vep -i ensembl-vep/inputs/allfam.vcf -o vep_anno.vcf --cache  --everything --dir_cache $HOME/cache "exact" --vcf

## bystro annotation

bystro annotation used website

## Filtering VEP annotation

VEP data is filtered and cleaned using the script vep_clean

## Filtering bystro Data

Script annotation.r is depreciated as project is pivoting to VEP;

