#script for cleaning and filtering bystro data
#used for CTCF family; script copied for newfam with new version bystro_annontation_newfam2.R
#===============================================================================
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse) #data cleaning
library(openxlsx) # excel
library(stringr) #string manipulation
library(vcfR) #vcf handling 
#READ in and first cleaning------------------------------------------------------
annotation_45 <- read.delim("~/Documents/Yang Lab/research/exsomes/allfam.annotation.tsv", header=TRUE,stringsAsFactors = F) 
#annotation_45<-read.delim("~/Documents/Yang Lab/research/exsomes/newfam2_vcf.annotation.tsv", header=TRUE,stringsAsFactors = F)

head_vcf <- read.table("~/Documents/Yang Lab/research/allfam.vcf", quote="\"") %>% 
  dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)

# head_vcf <- read.table("~/ensembl-vep/inputs/newfam2.vcf", quote="\"") %>% 
#   dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)

# Finder max values and determining variants where all transcripts are 
# intergenic with no function

max_cadd_finder<-str_split(annotation_45$cadd,"[;|]+")
annotation_45$max_cadd<-as.numeric(sapply(max_cadd_finder, max))

max_phyloP_finder<-str_split(annotation_45$phyloP,"[;|]+")
annotation_45$max_phyloP<-as.numeric(sapply(max_phyloP_finder, max))

max_phastCons_finder<-str_split(annotation_45$phastCons,"[;|]+")
annotation_45$max_phastCons<-as.numeric(sapply(max_phastCons_finder, max))

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
# filters based on allele frequency, quality, cadd, phylop, phastcons, allele function and site type
# additionally, add any entry with a non-benign entry in clinvar
pathogenic<-annotation_46 %>% filter(as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% 
  filter(max_cadd>15|max_phastCons>=0.9|max_phyloP>=3) %>% 
  filter(!(allele_func==TRUE&site_type==TRUE)) %>% 
  filter(!grepl('2002145',missingGenos)) %>% 
  filter(clinvar.check==TRUE) %>% 
  rbind(annotation_46 %>% filter(clinvar.check==FALSE,as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
          filter(clinvar.clinicalSignificance!="Benign"))


pathogenic$Mother<-ifelse(grepl('2002155',pathogenic$homozygotes),"Homozygote",
                             ifelse(grepl('2002155',pathogenic$heterozygotes),"Heterozygote","Not.Present"))
pathogenic$Father<-ifelse(grepl('2002161',pathogenic$homozygotes),"Homozygote",
                             ifelse(grepl('2002161',pathogenic$heterozygotes),"Heterozygote","Not.Present"))
pathogenic$Proband<-ifelse(grepl('2002145',pathogenic$homozygotes),"Homozygote","Heterozygote")

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

file_n<-paste(Sys.Date(),"bystro_family_2002.csv",sep = "_")

output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[+]","(+)")
output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[-]","(-)")
output_pathogenic$alt<-str_replace(output_pathogenic$alt,"[=]","")

write.csv(output_pathogenic,file = file_n,quote = F,row.names = F)





























