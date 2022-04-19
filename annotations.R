pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(tidyverse)
library(openxlsx)
library(stringr)
library(vcfR)
#READ in and first cleaning------------------------------------------------------
annotation_45 <- read.delim("~/Documents/Yang Lab/research/exsomes/allfam.annotation.tsv", header=TRUE,stringsAsFactors = F) 
# head_vcf<-read.vcfR("~/Documents/Yang Lab/research/allfam.vcf")
head_vcf <- read.table("~/Documents/Yang Lab/research/allfam.vcf", quote="\"") %>% 
  dplyr::select(chrom=V1,vcfPos=V2,ID=V3,QUAL=V6,FILTER=V7,FORMAT=V8)
#vcf_data_fields<-as.data.frame(head_vcf@fix) %>% rename(pos=POS,chrom=CHROM) %>% mutate(pos=as.numeric(pos))
#vcf_data_fields_2<-as.data.frame(head_vcf@gt)


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
pathogenic<-annotation_46 %>% filter(as.numeric(gnomad.genomes.af)<0.01|gnomad.genomes.af=="!") %>% 
  filter(as.numeric(QUAL)>20|is.na(QUAL)) %>% 
  #filter(max_cadd>30|max_phastCons>=0.9|max_phyloP>=3) %>%
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

allfam_vep <- read.table("~/ensembl-vep/outputs/allfam.tsv", header=T, quote="\"",stringsAsFactors = F) %>%
   dplyr::select(Location,cDNA_position) %>%
  separate(Location,into = c("chrom","chromStart"),sep = ":") %>% separate(chromStart,into=c("chromStart","chromEnd2"),sep="-") %>%
  mutate(chromEnd=as.numeric(chromStart)+1) %>% select(chrom,chromStart,cDNA_position)
allfam_vep2<-allfam_vep %>% filter(cDNA_position!="-") %>% group_by(chrom,chromStart) %>% 
  mutate(subrank=rank(cDNA_position,ties.method ="random")) %>% spread(value=cDNA_position,key = subrank,fill = "") %>%  unite(col = "cDNA",sep = ",") %>%
  separate(col="cDNA",into = c("chrom","chromStart","cDNA"),sep = ",",extra = "merge") %>% mutate(cDNA=str_remove_all(string = cDNA,pattern = ",,"))
allfam_vep2$chromStart2<-as.numeric(allfam_vep2$chromStart)+1
allfam_vep2$chromStart3<-as.numeric(allfam_vep2$chromStart)-1
allfam_vep2$chromStart4<-as.numeric(allfam_vep2$chromStart)+2
allfam_vep2$chromStart5<-as.numeric(allfam_vep2$chromStart)-2


allfam_vep2<-allfam_vep2 %>% gather(-chrom,-cDNA,key="buffer",value = "chromStart") %>% select(-buffer)

allfam_vep3<-allfam_vep2[(!duplicated(allfam_vep2[,c(1,3)])),]

output_pathogenic<-left_join(output_pathogenic,allfam_vep3,by=c("chrom","chromStart"))

output_pathogenic$cDNA<-ifelse(is.na(output_pathogenic$cDNA),"!",output_pathogenic$cDNA)

#Final ordering and output-------------------------------------------------------

col_order<-c("chrom","chromStart","chromEnd","ID","ref","alt","QUAL","cDNA","refSeq.refAminoAcid","refSeq.altAminoAcid","refSeq.codonNumber",
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
output_pathogenic<-output_pathogenic[,col_order]
output_pathogenic<-output_pathogenic %>% rename(clinvar.structure.phenotypeList=refSeq.clinvar.phenotypeList,
                             clinvar.structure.clinicalSignificance=refSeq.clinvar.clinicalSignificance,
                             clinvar.structure.type=refSeq.clinvar.type,                
                             clinvar.structure.origin=refSeq.clinvar.origin,
                             clinvar.structure.numberSubmitters=refSeq.clinvar.numberSubmitters,
                             clinvar.structure.reviewStatus=refSeq.clinvar.reviewStatus)



write.table(output_pathogenic,"b_annotation.bed",sep = "\t",quote = F,row.names = F)
 
#write_csv(sorted,"pathogenic_variant_file.csv")                                                                                                       
#----------------------------------------------------------------------------------






























