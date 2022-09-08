#Script for creating figure files for manuscript
# copies main script for filtering each source and creates figures
#==============================================================================
# Load needed libraries
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(openxlsx)
library(stringr)
library(vcfR)
library(ensemblVEP)
library(VariantAnnotation)


file<-"ensembl-vep/outputs/fa_allfam.vcf"
#file<-"ensembl-vep/outputs/newfam2_anon.vcf"
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

#=================================================================================
filtered_df<-merged_df%>%
  separate(PolyPhen,into = c("PolyPhen.1","PolyPhen.2"),sep = "[(]", remove = F) %>% 
  mutate(PolyPhen.2=str_sub(PolyPhen.2,start = 1,-2)) %>% dplyr::select(-PolyPhen.1) %>% 
  separate(SIFT,into = c("SIFT.1","SIFT.2"),sep = "[(]", remove = F) %>% 
  mutate(SIFT.2=str_sub(SIFT.2,start = 1,-2)) %>% dplyr::select(-SIFT.1) %>% 
  filter(as.numeric(gnomAD_AF)<0.001|is.na(gnomAD_AF),as.numeric(AF)<0.001|is.na(AF),as.numeric(MAX_AF)<0.05|is.na(MAX_AF)) %>%  
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% filter(as.numeric(SIFT.2<0.1)|is.na(SIFT)) %>% 
  filter(as.numeric(PolyPhen.2>0.9)|is.na(PolyPhen)) %>% 
  filter(!is.na(BIOTYPE),!is.na(Consequence)) %>%
  #filter(!is.na(`1871274`)) %>% 
  filter(!is.na(`2002145`)) %>%
  filter(!Consequence%in%c("downstream_gene_variant","upstream_gene_variant",
                           "non_coding_transcript_exon_variant","non_coding_transcript_variant","intergenic_variant",
                           "intron_variant","intron_variant&non_coding_transcript_variant","regulatory_region_variant",
                           "intron_variant&NMD_transcript_variant")) %>%  
  mutate_all(as.character)



fig_a<-filtered_df %>% group_by(IMPACT) %>%tally() %>% mutate(total=sum(n),Share=n/total,Type="Filtered") %>% 
  rbind(merged_df %>% group_by(IMPACT) %>%tally() %>% mutate(total=sum(n),Share=n/total,Type="Unfiltered")) %>% 
  ggplot(aes(x=Type,y=Share,fill=factor(IMPACT,levels = c("HIGH","MODERATE","LOW","MODIFIER"))))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_x_discrete(limits=c("Unfiltered","Filtered"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1),expand = c(0,0))+
  scale_fill_manual(breaks = c("HIGH","MODERATE","LOW","MODIFIER"),
                    values = c("#e31a1c","#fd8d3c","#fecc5c","#ffffb2"))+
  labs(title = "",x="",y="",
       fill="IMPACT")+theme_linedraw()+
  theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  coord_flip()
  

color_scale<-c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
              "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99")

fig_b<-filtered_df %>% group_by(BIOTYPE) %>%tally() %>% mutate(total=sum(n),Share=n/total,Type="Filtered") %>% 
  mutate(BIOTYPE.2=ifelse(base::grepl("pseudogene",BIOTYPE),"pseudogene",
                          ifelse(base::grepl("RNA",BIOTYPE),"micro or noncoding RNA",
                          ifelse(base::grepl("promoter|enhancer",BIOTYPE),"promoter/enhancer/promoter region",
                          ifelse(base::grepl("_gene",BIOTYPE),"immunoglobulin genes",
                          ifelse(base::grepl("decay|sense",BIOTYPE),"non-stop/nonsense-mediated decay/antisense",
                          ifelse(base::grepl("intron",BIOTYPE),"reatined intron/sense-intronic",
                          ifelse(base::grepl("binding_site",BIOTYPE),"CTCF/TF binding site",BIOTYPE)))))))) %>%  
  rbind(merged_df %>% group_by(BIOTYPE) %>%tally() %>% mutate(total=sum(n),Share=n/total,Type="Unfiltered") %>% 
          mutate(BIOTYPE.2=ifelse(base::grepl("pseudogene",BIOTYPE),"pseudogene",
                          ifelse(base::grepl("RNA",BIOTYPE),"micro or noncoding RNA",
                          ifelse(base::grepl("promoter|enhancer",BIOTYPE),"promoter/enhancer/promoter region",
                          ifelse(base::grepl("_gene",BIOTYPE),"immunoglobulin genes",
                          ifelse(base::grepl("decay|sense",BIOTYPE),"non-stop/nonsense-mediated decay/antisense",
                          ifelse(base::grepl("intron",BIOTYPE),"reatined intron/sense-intronic",
                          ifelse(base::grepl("binding_site",BIOTYPE),"CTCF/TF binding site",BIOTYPE))))))))) %>% 
  mutate(BIOTYPE.2=ifelse(is.na(BIOTYPE.2),"Other",BIOTYPE.2)) %>% 
  group_by(Type,BIOTYPE.2) %>% summarise(n=sum(n),Share=sum(Share)) %>% 
  mutate(BIOTYPE.2=factor(BIOTYPE.2,levels = c("protein_coding","immunoglobulin genes","CTCF/TF binding site","promoter/enhancer/promoter region","non-stop/nonsense-mediated decay/antisense",
                                               "pseudogene","micro or noncoding RNA","reatined intron/sense-intronic","open_chromatin_region","processed_transcript","Other"))) %>% 
  ggplot(aes(x=Type,y=Share,fill=BIOTYPE.2))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_x_discrete(limits=c("Unfiltered","Filtered"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1),expand = c(0,0))+
  scale_fill_manual(values = color_scale,
                   breaks = c("protein_coding","immunoglobulin genes","CTCF/TF binding site","promoter/enhancer/promoter region","non-stop/nonsense-mediated decay/antisense",
                              "pseudogene","micro or noncoding RNA","reatined intron/sense-intronic","open_chromatin_region","processed_transcript","Other"))+
  labs(title="",
       fill="BIOTYPE",x="",y="")+
  theme_linedraw()+
  theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  coord_flip()


cons<-filtered_df$Consequence %>% str_split("&") %>% unlist %>% base::as.data.frame() %>% mutate(Type="Filtered") %>% 
  rbind(merged_df$Consequence %>% str_split("&") %>% unlist %>% base::as.data.frame() %>% mutate(Type="Unfiltered"))

names(cons)<-c("Consequence","Type")

fig_c<-cons %>% group_by(Type,Consequence) %>% tally() %>% mutate(total=sum(n),Share=n/total) %>% 
  mutate(Consequence=ifelse(base::grepl("TFBS|TF_binding|NMD|miRNA|stream|UTR|elongation|regulatory",Consequence),"Other modiefier variants",
                            ifelse(base::grepl("inframe|frameshift|terminal",Consequence),"inframe insertion/deletion",
                            ifelse(base::grepl("non_coding",Consequence),"non-coding",
                            ifelse(base::grepl("splice",Consequence),"splice donor/acceptor",
                            ifelse(base::grepl("start|stop",Consequence),"loss-of-function",Consequence)))))) %>% 
  mutate(Consequence=factor(Consequence,
                            levels = c("loss-of-function","splice donor/acceptor","inframe insertion/deletion","missense_variant",
                                       "protein_altering_variant","coding_sequence_variant","synonymous_variant","non-coding","intron_variant",
                                       "intergenic_variant","Other modiefier variants"))) %>%
  ggplot(aes(x=Type,y=Share,fill=Consequence))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_x_discrete(limits=c("Unfiltered","Filtered"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1),expand = c(0,0))+
  scale_fill_manual(values = color_scale)+
  labs(title="",
       x="",
       y="")+
  theme_linedraw()+theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  coord_flip()




p_1<-filtered_df %>% dplyr::select(SIFT.2) %>% filter(!is.na(SIFT.2)) %>% mutate(Type="Filtered") %>% 
  rbind(merged_df %>% separate(SIFT,into = c("SIFT.1","SIFT.2"),sep = "[(]", remove = F) %>% 
          mutate(SIFT.2=str_sub(SIFT.2,start = 1,-2)) %>% dplyr::select(SIFT.2) %>% filter(!is.na(SIFT.2)) %>% mutate(Type="Unfiltered")) %>% mutate(SIFT.2=as.numeric(SIFT.2)) %>% 
  ggplot(aes(fill=Type,x=SIFT.2))+
  geom_density(alpha=0.5)+
  scale_y_continuous(limits = c(0,50))+
  scale_x_continuous(expand = c(0,0))+
  theme(legend.position = "none",
        title = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(title = "",x="SIFT")


p_2<-filtered_df %>% dplyr::select(PolyPhen.2) %>% filter(!is.na(PolyPhen.2)) %>% mutate(Type="Filtered") %>% 
  rbind(merged_df %>% separate(PolyPhen,into = c("PolyPhen.1","PolyPhen.2"),sep = "[(]", remove = F) %>% 
          mutate(PolyPhen.2=str_sub(PolyPhen.2,start = 1,-2)) %>% dplyr::select(PolyPhen.2) %>% filter(!is.na(PolyPhen.2)) %>% mutate(Type="Unfiltered")) %>% mutate(PolyPhen.2=as.numeric(PolyPhen.2)) %>% 
  ggplot(aes(fill=Type,x=PolyPhen.2))+
  geom_density(alpha=0.5)+
  scale_y_continuous(limits = c(0,50))+
  scale_x_continuous(expand = c(0,0))+
  theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(title = "",
      x="PolyPhen",y="")

#png(filename = "vep_scores.png",width = 800,height = 400)
cowplot::plot_grid(p_1,p_2,labels=c("",""))
ggsave(filename = "vep_scores.png",width = 8,height = 4,units = "in")
#dev.off()

#png(filename = "vep_image.png",width = 800,height = 1000)
cowplot::plot_grid(fig_a,fig_b,fig_c,nrow = 3,labels = "AUTO")
ggsave(filename = "vep_image.png",width = 8,height = 10,units = "in")
fig_a
#ggsave(filename = "fig_a.png",width = 8,height = 10,units = "in")
fig_b
#ggsave(filename = "fig_b.png",width = 8,height = 10,units = "in")
fig_c
#ggsave(filename = "fig_c.png",width = 8,height = 10,units = "in")
#cowplot::plot_grid(p_1,p_2,labels=c("",""))
#ggsave(filename = "fig_d.png",width = 8,height = 10,units = "in")

#dev.off()
rm(chrom_cross,cons,fig_a,fig_b,fig_c,fig_d,p_1,p_2)

#================================================================================
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)

annotation_45 <- read.delim("~/Documents/Yang Lab/research/exsomes/allfam.annotation.tsv", header=TRUE,stringsAsFactors = F) 
# head_vcf<-read.vcfR("~/Documents/Yang Lab/research/allfam.vcf")
head_vcf <- read.table("~/Documents/Yang Lab/research/allfam.vcf", quote="\"") %>% 
  dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)
#vcf_data_fields<-as.data.frame(head_vcf@fix) %>% rename(pos=POS,chrom=CHROM) %>% mutate(pos=as.numeric(pos))
#vcf_data_fields_2<-as.data.frame(head_vcf@gt)


max_cadd_finder<-str_split(annotation_45$cadd,"[;|]+")
suppressWarnings(annotation_45$max_cadd<-as.numeric(sapply(max_cadd_finder, max)))

max_phyloP_finder<-str_split(annotation_45$phyloP,"[;|]+")
suppressWarnings(annotation_45$max_phyloP<-as.numeric(sapply(max_phyloP_finder, max)))

max_phastCons_finder<-str_split(annotation_45$phastCons,"[;|]+")
suppressWarnings(annotation_45$max_phastCons<-as.numeric(sapply(max_phastCons_finder, max)))

# allele_func<-str_split(annotation_45$refSeq.exonicAlleleFunction, "[;|]+")
# annotation_45$allele_func<-sapply(allele_func, function(x){
#   (all(unlist(x)=="synonymous")|all(unlist(x)=="!"))
# }

allele_func<-str_split(annotation_45$refSeq.exonicAlleleFunction, "[;|]+")
annotation_45$allele_func<-sapply(allele_func, function(x){
  (all(unlist(x)=="!"))
})
site_type<-str_split(annotation_45$refSeq.siteType, "[;|]+")

annotation_45$site_type<-sapply(site_type, function(x){
  (all(unlist(x)=="intergenic"))
}) 

clinvar_check<-str_split(annotation_45$clinvar.alleleID,"[;|]+")
annotation_45$clinvar.check<-sapply(clinvar_check, function(x){
  (all(unlist(x)=="!"))
}) 
annotation_46<-left_join(annotation_45,head_vcf,by=c("chrom","vcfPos"))
#Filtering----------------------------------------------------------------------
suppressWarnings(pathogenic<-annotation_46 %>% filter(as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% 
  #filter(max_cadd>30|max_phastCons>=0.9|max_phyloP>=3) %>%
  filter(max_cadd>15|max_phastCons>=0.9|max_phyloP>=3) %>% 
  filter(!(allele_func==TRUE&site_type==TRUE)) %>% 
  filter(!grepl('2002145',missingGenos)) %>% 
  filter(clinvar.check==TRUE) %>% 
  rbind(annotation_46 %>% filter(clinvar.check==FALSE,as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
          filter(clinvar.clinicalSignificance!="Benign")))


funcs<-base::as.data.frame(unlist(allele_func)) %>% mutate(Type="Unfiltered")
funcs_2<-pathogenic$refSeq.exonicAlleleFunction %>% str_split("[;|]") %>% unlist() %>% base::as.data.frame() %>% mutate(Type="Filtered")
names(funcs)[1]<-"AlleleFunction"
names(funcs_2)[1]<-"AlleleFunction"

funcs_3<-annotation_45 %>% select(refSeq.exonicAlleleFunction,allele_func) %>% filter(allele_func==TRUE) %>% mutate(AlleleFunction="!",Type="Unfiltered") %>% select(AlleleFunction,Type)

color_scale_2<-c("#7fc97f","#fdc086","#ffff99","#386cb0","#bf5b17","#666666")

fig_e<-pathogenic %>% select(AlleleFunction=refSeq.exonicAlleleFunction) %>% mutate(Type="Filtered") %>% rbind(annotation_46 %>% 
                                                                                                                 select(AlleleFunction=refSeq.exonicAlleleFunction) %>% 
                                                                                                                 mutate(Type="Unfiltered")) %>% 
  mutate(AlleleFunction=ifelse(grepl("start|stop",AlleleFunction),"loss-of-function",
                               ifelse(grepl("nonSynonymous",AlleleFunction),"nonSynonymous",
                                      ifelse(grepl("indel-frameshift",AlleleFunction),"indel-frameshift",
                                             ifelse(grepl("indel-nonFrameshift",AlleleFunction),"indel-nonFrameshift",
                                                    ifelse(grepl("synonymous",AlleleFunction),"synonymous","no function")))))) %>% 
  group_by(Type,AlleleFunction) %>%  tally() %>% mutate(share=n/sum(n)) %>% 
  mutate(`Allele Function`=factor(AlleleFunction,
                                  levels = c("loss-of-function","nonSynonymous","indel-frameshift","indel-nonFrameshift","synonymous","no function"))) %>% 
  ggplot(aes(x=Type,fill=`Allele Function`,y=share))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_x_discrete(limits=c("Unfiltered","Filtered"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1),expand = c(0,0))+
  scale_fill_manual(values = color_scale_2)+
  theme_linedraw()+
  theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(x="",y="",title = "")+
  coord_flip()
  

site_1<-unlist(site_type) %>% as.data.frame() %>% mutate(Type="Unfiltered")
site_2<-unlist(str_split(pathogenic$refSeq.siteType,"[;|]")) %>% as.data.frame() %>% mutate(Type="Filtered")

names(site_1)[1]<-"Site"
names(site_2)[1]<-"Site"
color_scale_3<-c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17","#666666")

fig_f<-rbind(site_1,site_2) %>% group_by(Type,Site) %>% tally() %>%  mutate(share=n/sum(n)) %>% 
  mutate(Site=factor(Site,
                     levels = c("exonic","spliceDonor","spliceAcceptor","ncRNA","UTR3","UTR5","intronic","intergenic"))) %>% 
  ggplot(aes(x=Type,fill=Site,y=share))+
  geom_bar(stat = "identity",width = 0.5)+
  scale_x_discrete(limits=c("Unfiltered","Filtered"))+
  scale_y_continuous(labels = scales::percent,limits = c(0,1),expand = c(0,0))+
  scale_fill_manual(values = color_scale_3)+
  theme_linedraw()+
  theme(title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        axis.title=element_text(face = "bold"),
        axis.text = element_text(face = "bold"))+
  labs(x="",y="",title = "")+
  coord_flip()

#png(filename = "bystro_image.png",width = 800,height = 400)
cowplot::plot_grid(fig_f,fig_e,labels = c("A","B"),nrow = 2)
ggsave(filename = "bysto_image.png",width = 8,height = 6,units = "in")
#dev.off()

fig_e
#ggsave(filename = "fig_e.png",width = 8,height = 10,units = "in")
fig_f
#ggsave(filename = "fig_f.png",width = 8,height = 10,units = "in")

rm(head_vcf,max_cadd_finder,max_phastCons_finder,max_phyloP_finder,site_type)
rm(allele_func,clinvar_check,annotation_45)
 