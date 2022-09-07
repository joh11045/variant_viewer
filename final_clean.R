sorted <- read.delim("~/Documents/sorted.bed")  
pLI<-read.delim("~/Documents/gene_patho.bed")



sorted$oe<-ifelse(sorted$oe==".","!",sorted$oe)
sorted$pLI<-ifelse(sorted$pLI==".","!",sorted$pLI)  


median(as.numeric(sorted$phastCons),na.rm = T)
sd(as.numeric(sorted$phastCons),na.rm = T)

median(as.numeric(sorted$cadd),na.rm = T)
sd(as.numeric(sorted$cadd),na.rm = T)

median(as.numeric(sorted$phyloP),na.rm = T)
sd(as.numeric(sorted$phyloP),na.rm = T)

write_csv(sorted,file = "Variant_viewer/pathogenic_variants_report.csv")


allele_func2<-str_split(sorted$refSeq.exonicAlleleFunction, "[;|]+")
annotation_45$allele_func<-sapply(allele_func, function(x){
  (all(unlist(x)=="!"))
})

length(unlist(allele_func2))


site_1<-str_split(sorted$refSeq.siteType, "[;|]+")
site_2<-str_split(annotation_46$refSeq.siteType, "[;|]+")

length(unlist(site_1))
length(unlist(site_2))
