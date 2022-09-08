# Script for comparing output between VEP and bystro filtering
# not part of typlical workflow
library(tidyverse)

bystro<-read.table("2022-08-31_b_annotation.bed",header = T)

vep<-read.csv("family_18712.csv")

x<-sapply(bystro$ID, function(x){
  any(grepl(x,vep$Existing_variation))
})

y<-sapply(vep$Existing_variation, function(x){
  any(grepl(x,bystro$ID))
})


xx<-sapply(bystro$gene, function(x){
  any(grepl(x,vep$SYMBOL))
})

yy<-sapply(vep$SYMBOL, function(x){
  any(grepl(x,bystro$gene))
})
summary(bystro$ID%in%vep$Existing_variation)
summary(vep$Existing_variation%in%bystro$ID)


cadidates<-bystro %>%  
  group_by(gene) %>% tally() #%>% 
  left_join(bystro) %>% filter(Proband=="Heterozygote",sibling!="Not.Present",(Mother=="Not.Present"|Father=="Not.Present")) %>% filter(n>1)

#  write.csv(x,"8_10_candidates_bystro.csv",row.names = F)

x<-bystro %>% filter(Proband=="Heterozygote",sibling!="Not.Present",(Mother=="Not.Present"|Father=="Not.Present"),Mother!=Father) %>% 
  left_join(cadidates) %>% filter(n>1) %>%  group_by(gene) %>% 
  mutate(in_any_mom=all(Mother=="Not.Present"),in_any_dad=all(Father=="Not.Present")) %>% filter(in_any_mom==FALSE,in_any_dad==FALSE) %>% 
  dplyr::select(gene,n,Mother,Father,Proband,sibling,chrom,chromStart,chromEnd,ID,ref,alt,QUAL,refSeq.refAminoAcid,refSeq.altAminoAcid,refSeq.codonNumber,
             hgnc.gene.description,
             pubmed_links,omim_link,
             phastCons,phyloP,cadd,gnomad.genomes.af,gnomad_links,
             refSeq.nearest.name,
             refSeq.siteType,refSeq.exonicAlleleFunction,
             clinvar.alleleID,clinvar.clinicalSignificance,clinvar.type,
             clinvar.phenotypeList,clinvar.numberSubmitters,clinvar.origin, 
             clinvar.referenceAllele,clinvar.alternateAllele,clinvar.reviewStatus,        
             clinvar.structure.clinicalSignificance,clinvar.structure.type,clinvar.structure.phenotypeList,                
             clinvar.structure.numberSubmitters,clinvar.structure.origin,clinvar.structure.reviewStatus,refSeq.description,HGMD_link,genecard_links,db_best_guess)
x$alt<-str_replace(x$alt,"[+]","(+)")
x$alt<-str_replace(x$alt,"[-]","(-)")
x$alt<-str_replace(x$alt,"[=]","")

#---------------------------------------------------------------------------

#byst<-read.delim("~/Documents/variant_viewer/data/final_recode.annotation.tsv",header = T)
byst<-read.delim("~/Documents/Yang Lab/research/exsomes/newfam2_vcf.annotation.tsv",header=T)
file<-"ensembl-vep/outputs/newfam2_anon.vcf"
#file<-"ensembl-vep/outputs/final_newfam2.vcf"
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
                              ifelse(allele_1==0,"Heterozygous","Multiple variants")))) %>% 
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
filtered_df<-merged_df %>% mutate_all(as.character)

filtered_df[is.na(filtered_df)]<-"!"

#Simplify table to single row variant with multiple effects from each effect having
#its own row--------------------------------------------------------------------
frame_1<-filtered_df[!duplicated(filtered_df[,c("chr","pos","ref","alt")]),]
frame_2<-filtered_df[duplicated(filtered_df[,c("chr","pos","ref","alt")]),]

#individs<-c("2002145","2002155","2002161")
individs<-c("1871274","1871284","1871286","1871292")
for(i in 1:nrow(frame_1)){
  for(j in i:(i+95)){
    if(all(frame_1[i,c("chr","pos","ref","alt",individs,"Allele")]==frame_2[j,c("chr","pos","ref","alt",individs,"Allele")])){
      for(k in 9:ncol(frame_1)){
        frame_1[i,k]=ifelse(frame_1[i,k]==frame_2[j,k],frame_1[i,k],
                            base::paste(frame_1[i,k],frame_2[j,k],sep = "; "))
      }
    }
  }
}
library(dplyr)
x<-frame_1 %>% dplyr::select(chr,pos,Consequence,BIOTYPE) %>% mutate(pos=as.numeric(pos)) %>% 
  left_join(byst %>% dplyr::select(chr=chrom,pos,refSeq.siteType,refSeq.exonicAlleleFunction))

x$refSeq.siteType %>% table(useNA="ifany")

x$allele_func<-unlist(lapply(str_split(x$refSeq.exonicAlleleFunction, "[;|]+"),function(x){x[1]}))
x$Consequence_2<-unlist(lapply(str_split(x$Consequence,"[;&]"),function(x){x[1]}))
table(x$Consequence_2,x$allele_func,useNA="ifany")



a<-frame_1$Consequence %>% str_split(pattern = ';') %>% base::unlist() %>% str_split(pattern = "&") %>% unlist()
b<-byst$refSeq.siteType %>% str_split(pattern = ";") %>% base::unique()  %>%  base::unlist() %>% str_split(pattern = "[|]") %>% base::unique() %>% unlist()

table(a)
table(b)






output_pathogenic %>% inner_join(vep,by="key")


vep$key<-paste(vep$chr,vep$pos)
output_pathogenic$key<-paste(output_pathogenic$chr,output_pathogenic$pos)



vep<-X2022_09_07_vep_family_18712









merged_df$BIOTYPE %>% table()
merged_df$Consequence %>%str_split("[&]") %>% unlist() %>%  table()
annotation_46$refSeq.siteType %>% str_split("[|;]") %>% unlist() %>% table()
annotation_46$refSeq.exonicAlleleFunction %>% str_split("[|;]") %>% unlist() %>% table()


x<-unique(paste(merged_df$chr,merged_df$pos,sep = ":"))
y<-unique(paste(annotation_46$chrom,annotation_46$pos,sep = ":")) 

prop.table(table(x%in%y))
prop.table(table(y%in%x))

xx<-unique(paste(filtered_df$chr,filtered_df$pos,sep = ":"))
yy<-unique(paste(pathogenic$chrom,pathogenic$pos,sep = ":")) 

prop.table(table(xx%in%yy))
prop.table(table(yy%in%xx))

prop.table(table(xx%in%y))
prop.table(table(yy%in%x))

pathogenic$key<-paste(pathogenic$chrom,pathogenic$pos,sep = ":")
merged_df$key<-paste(merged_df$chr,merged_df$pos,sep = ":")

annotation_46$key<-paste(annotation_46$chrom,annotation_46$pos,sep = ":")

zz<-annotation_46 %>% filter(allele_func==T,site_type==T) %>% select(key,refSeq.siteType,refSeq.exonicAlleleFunction,ref,alt) %>% inner_join(merged_df %>% select(key,BIOTYPE,Consequence,ref,alt))

zzz<-zz %>% group_by(Consequence,BIOTYPE) %>% tally()

filtered_df$key<-paste(filtered_df$chr,filtered_df$pos,sep = ":")

z<-pathogenic %>% select(key,refSeq.siteType,refSeq.exonicAlleleFunction,ref,alt) %>% inner_join(merged_df %>% select(key,BIOTYPE,Consequence,ref,alt))

zz<-merged_df %>% filter(key%in%yy) %>% anti_join(filtered_df,by="key")
xxx<-pathogenic %>% filter(key%in%x) %>% anti_join(filtered_df,by="key")

