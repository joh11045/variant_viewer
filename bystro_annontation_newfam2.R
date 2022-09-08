#script for cleaning and filtering bystro data
#used for new family; script copied from annotations.R 
#===============================================================================
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(openxlsx)
library(stringr)
library(vcfR)
#READ in and first cleaning------------------------------------------------------
#annotation_45 <- read.delim("~/Documents/Yang Lab/research/exsomes/allfam.annotation.tsv", header=TRUE,stringsAsFactors = F) 
annotation_45<-read.delim("~/Documents/Yang Lab/research/exsomes/newfam2_vcf.annotation.tsv", header=TRUE,stringsAsFactors = F)
# head_vcf<-read.vcfR("~/Documents/Yang Lab/research/allfam.vcf")
# head_vcf <- read.table("~/Documents/Yang Lab/research/allfam.vcf", quote="\"") %>% 
#  dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)

head_vcf <- read.table("~/ensembl-vep/inputs/newfam2.vcf", quote="\"") %>% 
  dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)


max_cadd_finder<-str_split(annotation_45$cadd,"[;|]+")
annotation_45$max_cadd<-as.numeric(sapply(max_cadd_finder, max))

max_phyloP_finder<-str_split(annotation_45$phyloP,"[;|]+")
annotation_45$max_phyloP<-as.numeric(sapply(max_phyloP_finder, max))

max_phastCons_finder<-str_split(annotation_45$phastCons,"[;|]+")
annotation_45$max_phastCons<-as.numeric(sapply(max_phastCons_finder, max))

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
pathogenic<-annotation_46 %>% filter(as.numeric(gnomad.genomes.af)<0.001|gnomad.genomes.af=="!") %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% 
  #filter(max_cadd>30|max_phastCons>=0.9|max_phyloP>=3) %>%
  filter(max_cadd>15|max_phastCons>=0.9|max_phyloP>=3) %>% 
  filter(!(allele_func==TRUE&site_type==TRUE)) %>% 
  filter(!grepl('1871274',missingGenos)) %>% 
  filter(clinvar.check==TRUE) %>% 
  rbind(annotation_46 %>% filter(clinvar.check==FALSE,as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
          filter(clinvar.clinicalSignificance!="Benign"))


pathogenic$Mother<-ifelse(grepl('1871284',pathogenic$homozygotes),"Homozygote",
                          ifelse(grepl('1871284',pathogenic$heterozygotes),"Heterozygote","Not.Present"))
pathogenic$Father<-ifelse(grepl('1871286',pathogenic$homozygotes),"Homozygote",
                          ifelse(grepl('1871286',pathogenic$heterozygotes),"Heterozygote","Not.Present"))
pathogenic$Proband<-ifelse(grepl('1871274',pathogenic$homozygotes),"Homozygote","Heterozygote")

pathogenic$sibling<-ifelse(grepl('1871292',pathogenic$homozygotes),"Homozygote",
                           ifelse(grepl('1871292',pathogenic$heterozygotes),"Heterozygote","Not.Present"))
#hgnc_gene descriptions added---------------------------------------------------

hgnc_complete_set <- read.delim("~/Documents/Yang Lab/research/hgnc_complete_set.txt",stringsAsFactors = F)

ref_accession<-str_split(pathogenic$refSeq.name,"[;|]+")

w<-vector(mode = "list")
for (i in 1:length(pathogenic$chrom)) {
  w[[i]]<-unique(unlist(simplify2array(lapply(ref_accession[[i]], function(x){
    grep(pattern = x,hgnc_complete_set$refseq_accession)
  }))))}

gene_desc<-lapply(w,function(x){
  p<-list()
  for (j in 1:length(x)){
    p[j]<-hgnc_complete_set$name[x[j]]
  }
  sapply(p, paste)
})

gene_desc<-sapply(gene_desc, function(x){
  paste(unlist(x),collapse = "|")
})

pathogenic$hgnc.gene.description<-ifelse(gene_desc=="NA","!",gene_desc)


#cleaning/sorting---------------------------------------------------------------

output_pathogenic<-pathogenic
output_pathogenic$chrom_sort<-str_sub(output_pathogenic$chrom,start = 4)
output_pathogenic$chrom_sort<-ifelse(output_pathogenic$chrom_sort=="X",24,
                                     ifelse(output_pathogenic$chrom_sort=="Y",25,output_pathogenic$chrom_sort))
output_pathogenic<-output_pathogenic[order(output_pathogenic$chrom_sort,output_pathogenic$pos),]


output_pathogenic$cadd<-output_pathogenic$max_cadd
output_pathogenic$phastCons<-output_pathogenic$max_phastCons
output_pathogenic$phyloP<-output_pathogenic$max_phyloP

output_pathogenic<-output_pathogenic %>% mutate_at(c("cadd","phyloP","phastCons","gnomad.genomes.af","gnomad.exomes.af"),as.numeric)
output_pathogenic<-output_pathogenic %>% mutate_at(c("cadd","phyloP","phastCons","gnomad.genomes.af","gnomad.exomes.af"), ~round(.,digits = 5))
gene_names<-sapply(str_split(output_pathogenic$refSeq.name2,"[;|]+"),function(x){
  x[1]
})
db_snp<-sapply(str_split(output_pathogenic$dbSNP.name,"[;|]+"),function(x){
  x[1]
})
output_pathogenic$gene<-gene_names
output_pathogenic$ID<-db_snp
output_pathogenic$chromStart<-output_pathogenic$pos
output_pathogenic$chromEnd<-output_pathogenic$pos+1

#OMIM data----------------------------------------------------------------------
OMIM<-read.delim("~/Documents/Yang Lab/research/mim2gene.txt", header=FALSE, comment.char="#", stringsAsFactors=FALSE) 
names(OMIM)<-c("MIM_Number","MIM_Entry Type", "Entrez_Gene_ID", "Approved_Gene_Symbol", "Ensemble_Gene_ID")

OMIM<-OMIM %>% select(MIM_Number,Approved_Gene_Symbol) %>% filter(!duplicated(Approved_Gene_Symbol),Approved_Gene_Symbol!="")

output_pathogenic<-left_join(output_pathogenic,OMIM,by=c("gene"="Approved_Gene_Symbol")) %>% filter()




#output_pathogenic %>% mutate(chromStart=) %>% select(chrom)
#output_pathogenic %>% mutate(chrom=substring(chrom,4)) %>% dplyr::select(chrom,pos,ID) %>% 
# write.table("output_pathogenic.txt",row.names = F)

MIM<-output_pathogenic$MIM_Number

#Prep links--------------------------------------------------------------------
output_pathogenic$genecard_links<-ifelse(is.na(gene_names),"!",paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene_names))
output_pathogenic$pubmed_links<-ifelse(is.na(gene_names),"!",paste0("https://pubmed.ncbi.nlm.nih.gov/?term=",gene_names))
output_pathogenic$omim_link<-ifelse(is.na(MIM),"!",paste0("https://www.omim.org/entry/",MIM))
output_pathogenic$db_best_guess<-ifelse(is.na(db_snp)|db_snp=="!","!",paste0("https://www.ncbi.nlm.nih.gov/snp/",db_snp))
output_pathogenic$gnomad_links<-ifelse(is.na(db_snp),"!",paste0("https://gnomad.broadinstitute.org/variant/",db_snp))
output_pathogenic$HGMD_link<-ifelse(is.na(gene_names),"!",paste0("http://www.hgmd.cf.ac.uk/ac/gene.php?gene=",gene_names))

output_pathogenic<-as.data.frame(sapply(names(output_pathogenic),function(x){
  output_pathogenic[,x]<-ifelse(is.na(output_pathogenic[,x]),"!",output_pathogenic[,x])
})) 

output_pathogenic<-as.data.frame(sapply(names(output_pathogenic),function(y){
  sapply(str_split(output_pathogenic[,y],pattern =  ";|\\|"),function(x){
    str_c(unique(unlist(x)),collapse = ";")
  })}))

output_pathogenic<-as.data.frame(sapply(names(output_pathogenic),function(y){
  sapply(str_split(output_pathogenic[,y],pattern =  " "),function(x){
    str_c(unique(unlist(x)),collapse = "-")
  })}))

output_pathogenic$alt<-ifelse(grepl("=",output_pathogenic$alt),as.character(paste0(":",output_pathogenic$alt)),
                              as.character(output_pathogenic$alt))
#CDNA-------------------------------------------------------------------------------


#Final ordering and output-------------------------------------------------------
output_pathogenic$chr=output_pathogenic$chrom
output_pathogenic$pos=output_pathogenic$chromStart
col_order<-c("chr","pos","ref","alt","QUAL","gene","Proband","Mother","Father","sibling",
             "refSeq.refAminoAcid","refSeq.altAminoAcid","refSeq.codonNumber",
             "phastCons","phyloP","cadd","gnomad.genomes.af",
             "refSeq.nearest.name",
             "refSeq.siteType","refSeq.exonicAlleleFunction", "hgnc.gene.description",
             "clinvar.alleleID","clinvar.clinicalSignificance","clinvar.type",
             "clinvar.phenotypeList","clinvar.numberSubmitters","clinvar.origin", 
             "clinvar.referenceAllele","clinvar.alternateAllele","clinvar.reviewStatus",        
             "refSeq.clinvar.clinicalSignificance","refSeq.clinvar.type","refSeq.clinvar.phenotypeList",                
             "refSeq.clinvar.numberSubmitters","refSeq.clinvar.origin","refSeq.clinvar.reviewStatus","refSeq.description","gnomad_links",
             "pubmed_links","omim_link","HGMD_link","genecard_links","db_best_guess")
output_pathogenic<-output_pathogenic[,col_order]
output_pathogenic<-output_pathogenic %>% rename(clinvar.structure.phenotypeList=refSeq.clinvar.phenotypeList,
                                                clinvar.structure.clinicalSignificance=refSeq.clinvar.clinicalSignificance,
                                                clinvar.structure.type=refSeq.clinvar.type,                
                                                clinvar.structure.origin=refSeq.clinvar.origin,
                                                clinvar.structure.numberSubmitters=refSeq.clinvar.numberSubmitters,
                                                clinvar.structure.reviewStatus=refSeq.clinvar.reviewStatus)

file_n<-paste(Sys.Date(),"bystro_family_18712.csv",sep = "_")

output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[+]","(+)")
output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[-]","(-)")
output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[=]","")

write.csv(output_pathogenic,file = file_n,quote = F,row.names = F)

#write_csv(sorted,"pathogenic_variant_file.csv")                                                                                                       
#----------------------------------------------------------------------------------




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






























