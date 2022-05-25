# Load needed libraries
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(openxlsx)
library(stringr)
library(vcfR)
library(ensemblVEP)
library(VariantAnnotation)
#===============================================================================
#Clean up annotation data-------------------------------------------------------

#Read in data-------------------------------------------------------------------
## set file to VEP output, current process uses output from the VEP command
## ensembl-vep/./vep -i INPUTFILE -o OUTPUTFILE --cache --vcf --everything --dir_cache CACHE DIR --offline

#file<-"~/vep_anno.vcf"

file<-"ensembl-vep/outputs/allfamhg37.vcf"

## use functions from Variant Annotation Package to load VCF data in order to use
## functions from ensemblVEP package
trial<-readVcf(file)
trial.2<-parseCSQToGRanges(trial)

## CSQ refers to ConSeQuence and is the annotated data

csq<-as.data.frame(cbind(trial.2@ranges@NAMES,trial.2@elementMetadata)) %>%  
  separate(trial.2.ranges.NAMES,into = c("chr","temp"),sep = "[:]",extra = "merge") %>% 
  separate(temp,into = c("pos","ref_alt"),sep = "[_]")

rm(trial,trial.2)

#saveRDS(csq,"temp.rds")
#csq<-readRDS("temp.rds")

## Read in the same VCF file again but to identify actual alleles at each location and QUAL
allfam<-read.vcfR(file)
x<-vcfR2tidy(allfam,toss_INFO_column = T)

head_vcf<-x$fix
samps<-x$gt

## ChromKey does not correspond to chrom number, just sub in chr and avoid
chrom_cross<-head_vcf %>% group_by(ChromKey,chr=CHROM) %>% tally() %>% dplyr::select(-n)
samps<-left_join(samps,chrom_cross)

rm(x,allfam)
#Determine hetero/homozygosity--------------------------------------------------
## variable gt_GT list 0 for ref, 1 for alt, 1+n for each additional alt at location
## if allele is na; variant not present in individual
## if allele 1 = allele 2; homozygous. no individuals are 0/0, 
## if allele 1 = 0, then heterozygous, 
## if allele 1 dne 0 and dne allele 2 (rare); than individual has two variants at location
## current data is each individual per row, transform to one variant multiple individuals per row
de_novo_finder<-samps %>% dplyr::select(chr,POS,Indiv,gt_GQ,gt_DP,gt_GT,gt_GT_alleles) %>% 
  separate(gt_GT,into = c("allele_1","allele_2"),sep = "/") %>%
  mutate(status=ifelse(is.na(allele_1),"No Variant",
                       ifelse(allele_1==allele_2,"Homozygous",
                              ifelse(allele_1==0,"Heterozygous","Homozygous multiple variants")))) %>% 
  dplyr::select(chr,POS,gt_GT_alleles,Indiv,status) %>%
  filter(gt_GT_alleles!=".") %>% 
  base::unique() %>%
  spread(Indiv,status) 
 
#isolate quality scores---------------------------------------------------------
quals<-head_vcf %>% dplyr::select(chr=CHROM,pos=POS,ref=REF,alt=ALT,QUAL)

#merge annotations, vcf header, and samples-------------------------------------
temp_df<-de_novo_finder %>% dplyr::rename(pos=POS) %>% 
  mutate(pos=as.character(pos)) %>% 
  left_join(csq,by=c("chr","pos")) %>% separate(ref_alt,c("ref","alt"),sep='[/]')

merged_df<-temp_df %>% mutate(pos=as.numeric(pos)) %>% left_join(quals,by = c("chr","pos","ref","alt"))

rm(head_vcf,csq,de_novo_finder,quals,samps,temp_df)
#filter-------------------------------------------------------------------------
## isolate numbers from polyphen and sift columns
## filter based on SIFT polyPhen, allele frequency, QUALITY, BIOTYPE and Consequence
## Note: some variants usually complicated insertions or deletions will have no associated quality

filtered_df<-merged_df%>%
  separate(PolyPhen,into = c("PolyPhen.1","PolyPhen.2"),sep = "[(]", remove = F) %>% 
  mutate(PolyPhen.2=str_sub(PolyPhen.2,start = 1,-2)) %>% dplyr::select(-PolyPhen.1) %>% 
  separate(SIFT,into = c("SIFT.1","SIFT.2"),sep = "[(]", remove = F) %>% 
  mutate(SIFT.2=str_sub(SIFT.2,start = 1,-2)) %>% dplyr::select(-SIFT.1) %>% 
  filter(as.numeric(AF<0.01)|is.na(AF)) %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% filter(as.numeric(SIFT.2<0.1)|is.na(SIFT)) %>% 
  filter(as.numeric(PolyPhen.2>0.9)|is.na(PolyPhen)) %>% 
  filter(!is.na(BIOTYPE),!is.na(Consequence)) %>%
  filter(!is.na(`2002145`)) %>% 
  filter(IMPACT=="HIGH"|IMPACT=="MODERATE") %>% 
  mutate_all(as.character)

filtered_df[is.na(filtered_df)]<-"!"

#Simplify table to single row variant with multiple effects from each effect having
#its own row--------------------------------------------------------------------
frame_1<-filtered_df[!duplicated(filtered_df[,c("chr","pos","ref","alt")]),]
frame_2<-filtered_df[duplicated(filtered_df[,c("chr","pos","ref","alt")]),]

for(i in 1:nrow(frame_1)){
  for(j in i:(i+95)){
    if(all(frame_1[i,c("chr","pos","ref","alt","2002145","2002155","2002161","Allele")]==frame_2[j,c("chr","pos","ref","alt","2002145","2002155","2002161","Allele")])){
      for(k in 9:ncol(frame_1)){
        frame_1[i,k]=ifelse(frame_1[i,k]==frame_2[j,k],frame_1[i,k],
                            base::paste(frame_1[i,k],frame_2[j,k],sep = "; "))
    }
    }
  }
}
rm(chrom_cross,filtered_df,frame_2,merged_df,i,j,k)

col_order<-c("chr","ref","ID","ref","alt","QUAL","cDNA","refSeq.refAminoAcid","refSeq.altAminoAcid","refSeq.codonNumber",
             "gene","hgnc.gene.description","Mother","Father","Proband",
             "pubmed_links","omim_link",
             "phastCons","phyloP","cadd","gnomad.genomes.af","gnomad_links",
             "refSeq.nearest.name",
             "refSeq.siteType","refSeq.exonicAlleleFunction",
             "clinvar.alleleID","clinvar.clinicalSignificance","clinvar.type",
             "clinvar.phenotypeList","clinvar.numberSubmitters","clinvar.origin", 
             "clinvar.referenceAllele","clinvar.alternateAllele","clinvar.reviewStatus",        
             "refSeq.clinvar.clinicalSignificance","refSeq.clinvar.type","refSeq.clinvar.phenotypeList",                
             "refSeq.clinvar.numberSubmitters","refSeq.clinvar.origin","refSeq.clinvar.reviewStatus","refSeq.description","HGMD_link","genecard_links","db_best_guess")
