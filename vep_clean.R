pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(openxlsx)
library(stringr)
library(vcfR)
library(ensemblVEP)
library(VariantAnnotation)
#===============================================================================
#Read in data-------------------------------------------------------------------
#Clean up annotation data-------------------------------------------------------
trial<-readVcf("~/vep_anno.vcf",genome = "hg19")
trial.2<-parseCSQToGRanges(trial)


csq<-as.data.frame(cbind(trial.2@ranges@NAMES,trial.2@elementMetadata)) %>%  
  separate(trial.2.ranges.NAMES,into = c("chr","temp"),sep = "[:]",extra = "merge") %>% 
  separate(temp,into = c("pos","ref_alt"),sep = "[_]")

rm(trial,trial.2)

saveRDS(csq,"temp.rds")
csq<-readRDS("temp.rds")

allfam<-read.vcfR("~/vep_anno.vcf")


#samples<-(allfam@gt)
#allfam_vcf_annon<-(allfam@fix)

x<-vcfR2tidy(allfam,toss_INFO_column = T)

annotations<-x$fix
samps<-x$gt
rm(x)
#Determine hetero/homozygosity--------------------------------------------------
de_novo_finder<-samps %>% dplyr::select(ChromKey,POS,Indiv,gt_GQ,gt_DP,gt_GT,gt_GT_alleles) %>% separate(gt_GT,into = c("allele_1","allele_2"),sep = "/") %>%
  mutate(status=ifelse(is.na(allele_1),"No Variant",
                       ifelse(allele_1==allele_2,"Homozygous",
                              ifelse(allele_1==0,"Heterozygous","Homozygous multiple variants")))) %>% 
  dplyr::select(ChromKey,POS,gt_GT_alleles,Indiv,status) %>%
  filter(gt_GT_alleles!=".") %>% 
  base::unique() %>%
  spread(Indiv,status) 
 
de_novo_finder %>% head()
csq %>% head()

quals<-annotations %>% dplyr::select(chr=CHROM,pos=POS,ref=REF,alt=ALT,QUAL)
#merge annotations, vcf header, and samples-------------------------------------
temp_df<-de_novo_finder %>%  dplyr::rename(chr=ChromKey,pos=POS,ref_alt=gt_GT_alleles) %>%
  mutate(chr=paste0("chr",chr),pos=as.character(pos)) %>% left_join(csq) %>% separate(ref_alt,c("ref","alt"),sep='[/]')

temp_df %>% mutate(pos=as.numeric(pos)) %>% left_join(quals)->merged_df

#filter-------------------------------------------------------------------------
merged_df %>% filter(as.numeric(gnomAD_AF<0.01)|is.na(gnomAD_AF)) %>% 
  filter(QUAL>20) %>% 
  


pathogenic<-annotation_46 %>% filter(as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% 
  #filter(max_cadd>30|max_phastCons>=0.9|max_phyloP>=3) %>%
  filter(max_cadd>15|max_phastCons>=0.9|max_phyloP>=3) %>% 
  filter(!(allele_func==TRUE&site_type==TRUE)) %>% 
  filter(!grepl('2002145',missingGenos)) %>% 
  filter(clinvar.check==TRUE) %>% 
  rbind(annotation_46 %>% filter(clinvar.check==FALSE,as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
          filter(clinvar.clinicalSignificance!="Benign"))