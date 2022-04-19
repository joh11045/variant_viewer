# variant_viewer
Code for the development of a pathogenic variant viewer app
Variant viewer app displays table of possibly pathogenic variants and an IGV instance

Table is obtained through annotation and filtering. Filtering scheme is present in annotations.r

## VEP annotation

annotation was command line

ensembl-vep/./vep -i ensembl-vep/inputs/allfam.vcf -o vep_anno.vcf --cache  --everything --dir_cache $HOME/cache "exact" --vcf

## bystro annotation

bystro annotation used website
