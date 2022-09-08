## APP
Main app is variant_viewer_app.Rmd
includes newer feature allowing to filter by family variant status. Older app also included without that functionality

mardown_app.Rmd

both apps need files to be set to work properly

## VEP annotation

annotation was command line


ensembl-vep/./vep -i ensembl-vep/inputs/allfam.vcf -o vep_anno.tsv --cache  --everything --dir_cache $HOME/cache "exact" --tab

ensembl-vep/./vep -i ensembl-vep/inputs/allfam.vcf -o vep_anno.vcf --cache  --everything --dir_cache $HOME/cache "exact" --vcf

## bystro annotation

bystro annotation used website

## Filtering VEP annotation

VEP data is filtered and cleaned using the script vep_clean

## Filtering bystro Data

Script annotation.r or bystro_newfam2_annotation

