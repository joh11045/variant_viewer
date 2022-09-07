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

#file<-"ensembl-vep/outputs/fa_allfam.vcf"
#file<-"ensembl-vep/outputs/final_allfam.vcf"
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

filtered_df<-merged_df%>%
  separate(PolyPhen,into = c("PolyPhen.1","PolyPhen.2"),sep = "[(]", remove = F) %>% 
  mutate(PolyPhen.2=str_sub(PolyPhen.2,start = 1,-2)) %>% dplyr::select(-PolyPhen.1) %>% 
  separate(SIFT,into = c("SIFT.1","SIFT.2"),sep = "[(]", remove = F) %>% 
  mutate(SIFT.2=str_sub(SIFT.2,start = 1,-2)) %>% dplyr::select(-SIFT.1) %>% 
  filter(as.numeric(gnomAD_AF)<0.001|is.na(gnomAD_AF),as.numeric(AF)<0.001|is.na(AF),as.numeric(MAX_AF)<0.05|is.na(MAX_AF)) %>%  
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% filter(as.numeric(SIFT.2<0.1)|is.na(SIFT)) %>% 
  filter(as.numeric(PolyPhen.2>0.9)|is.na(PolyPhen)) %>% 
  filter(!is.na(BIOTYPE),!is.na(Consequence)) %>%
  filter(!is.na(`1871274`)) %>% 
  #filter(!is.na(`2002145`)) %>%
  filter(!Consequence%in%c("downstream_gene_variant","upstream_gene_variant",
                           "non_coding_transcript_exon_variant","non_coding_transcript_variant","intergenic_variant",
                           "intron_variant","intron_variant&non_coding_transcript_variant","regulatory_region_variant",
                           "intron_variant&NMD_transcript_variant")) %>%
  mutate_all(as.character)

# merged_df%>%separate(PolyPhen,into = c("PolyPhen.1","PolyPhen.2"),sep = "[(]", remove = F) %>% 
#   mutate(PolyPhen.2=str_sub(PolyPhen.2,start = 1,-2)) %>% dplyr::select(-PolyPhen.1) %>% 
#   separate(SIFT,into = c("SIFT.1","SIFT.2"),sep = "[(]", remove = F) %>% 
#   mutate(SIFT.2=str_sub(SIFT.2,start = 1,-2)) %>% dplyr::select(-SIFT.1) %>% 
#   filter(as.numeric(AF<0.01)|is.na(AF)) %>% 
#   filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% filter(as.numeric(SIFT.2<0.1)|is.na(SIFT)) %>% 
#   filter(as.numeric(PolyPhen.2>0.9)|is.na(PolyPhen)) %>% 
#   filter(!is.na(BIOTYPE),!is.na(Consequence)) %>%
#   filter(!is.na(`2002145`)) %>% 
#   filter(IMPACT!="LOW") %>%
#   filter(Consequence!="downstream_gene_variant",Consequence!="upstream_gene_variant",
#          Consequence!="non_coding_transcript_exon_variant",Consequence!="non_coding_transcript_variant",
#          Consequence!="intron_variant&non_coding_transcript_variant",Consequence!="intron_variant",
#          Consequence!="3_prime_UTR_variant",Consequence!="5_prime_UTR_variant",
#          Consequence!="intron_variant&NMD_transcript_variant")->y
# %>%
# merged_df %>%   group_by(Consequence,BIOTYPE,IMPACT) %>% tally()->x

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
frame_1<-frame_1 %>% filter(!(`1871274`=="Homozygous"&`1871284`=="!"&`1871286`=="!"&`1871292`=="!"))
rm(chrom_cross,filtered_df,frame_2,merged_df,i,j,k)
# Create links------------------------------------------------------------------
db_snp<-word(frame_1$Existing_variation,1,sep="&")
OMIM<-read.delim("~/Documents/Yang Lab/research/mim2gene.txt", header=FALSE, comment.char="#", stringsAsFactors=FALSE) 
names(OMIM)<-c("MIM_Number","MIM_Entry Type", "Entrez_Gene_ID", "Approved_Gene_Symbol", "Ensemble_Gene_ID")

OMIM<-OMIM %>% dplyr::select(MIM_Number,Approved_Gene_Symbol) %>% filter(!duplicated(Approved_Gene_Symbol),Approved_Gene_Symbol!="")

frame_1<-left_join(frame_1,OMIM,by=c("SYMBOL"="Approved_Gene_Symbol"))


frame_1$genecard_links<-ifelse(is.na(frame_1$SYMBOL),"!",paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",frame_1$SYMBOL))
frame_1$pubmed_links<-ifelse(is.na(frame_1$SYMBOL),"!",paste0("https://pubmed.ncbi.nlm.nih.gov/?term=",frame_1$SYMBOL))
frame_1$omim_link<-ifelse(is.na(frame_1$MIM_Number),"!",paste0("https://www.omim.org/entry/",frame_1$MIM_Number))
frame_1$db_best_guess<-ifelse(is.na(db_snp)|db_snp=="!","!",paste0("https://www.ncbi.nlm.nih.gov/snp/",db_snp))
frame_1$gnomad_links<-ifelse(is.na(db_snp),"!",paste0("https://gnomad.broadinstitute.org/variant/",db_snp))
frame_1$HGMD_link<-ifelse(is.na(frame_1$SYMBOL),"!",paste0("http://www.hgmd.cf.ac.uk/ac/gene.php?gene=",frame_1$SYMBOL))



final_table<-frame_1 %>% dplyr::select(chr, pos, ref, alt, QUAL,gene=SYMBOL,Proband=`1871274`,Mother=`1871284`,Father=`1871286`,Sibling=`1871292`,   #Proband=`2002145`, Mother=`2002155`, Father=`2002161`,
                                       Amino_acids, Codons,SIFT, PolyPhen, AF, gnomAD_AF, MAX_AF, MAX_AF_POPS,
                                       `HGVS Simple`=HGVSc, `HGVS Protien`=HGVSp,
                                       Consequence, BIOTYPE, IMPACT, Feature, Feature_type, EXON, INTRON, HGVSc, HGVSp, cDNA_position,
                                       CDS_position, Protein_position,  Existing_variation, DISTANCE, STRAND,     
                                       FLAGS, VARIANT_CLASS, SYMBOL_SOURCE, HGNC_ID, CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL,  
                                       TSL, APPRIS, CCDS, ENSP, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, GENE_PHENO,
                                        DOMAINS, miRNA, CLIN_SIG,             
                   SOMATIC, PHENO, PUBMED,HGMD_link,gnomad_links,db_best_guess,omim_link,pubmed_links,genecard_links)

final_table$alt<-str_replace(final_table$alt,"[+]","(+)")
final_table$alt<-str_replace(final_table$alt,"[-]","(-)")
final_table$alt<-str_replace(final_table$alt,"[=]","")

#denovos<-final_table %>% filter((Proband=="Homozygous"&Sibling=="Heterozygous")|(Proband=="Heterozygous"&Sibling=="!"))
#saveRDS(final_table,"Variant_viewer/vep_table.rds")

#x<-read_rds("Variant_viewer/vep_table.rds")
#x[1,] %>% gather(value,variable) %>% write.csv("Variant_viewer/vep_data_dict.csv",row.names = F)

#
#frame_1 %>% dplyr::select(chr,pos) -> x

#write_tsv(x,file = "ft_newfam_positions.txt")
file_n<-paste(Sys.Date(),"vep_family_18712.csv",sep = "_")

write.csv(final_table,file = file_n,quote = F,row.names = F)
