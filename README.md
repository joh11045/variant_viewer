# variant_viewer
Code for the development of a pathogenic variant viewer app
Variant viewer app displays table of possibly pathogenic variants and an IGV instance

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

