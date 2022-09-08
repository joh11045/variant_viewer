#Script to find variants for given pattern: multiple variants on a gene and with a pattern of missing in one parent and present in the other for each parent 
final_table<-read.csv("family_18712.csv")
denovos<-read.csv("family_18712_not_in_sibling.csv")

x<-base::as.data.frame(unlist(stringr::str_split(string = frame_1$Consequence,pattern  = "&")))
names(x)<-"Consequence"
x<-x %>% group_by(Consequence) %>% tally() %>% arrange(n)
ggplot(x)+
  geom_bar(aes(x=reorder(Consequence, n, sum),y=n),stat = "identity")+
  labs(x="Consequence",title = "Filtered Variants by Consequence")+
  coord_flip()

  summary.factor(x)

summary.factor(frame_1$IMPACT)
summary.factor(denovos$IMPACT)

summary.factor(frame_1$BIOTYPE)

ggplot(frame_1)+
  geom_histogram(aes(x=as.numeric(SIFT.2)))

frame_1 %>% mutate(SIFT.3=ifelse(SIFT.2=="!","Missing",ifelse(SIFT.2<0.05,"Deleterious","Tolerated"))) %>% group_by(SIFT.3) %>% tally()

subset<-final_table %>% group_by(SYMBOL) %>% summarise(variants_on_gene=n()) %>% left_join(final_table) %>% filter(variants_on_gene>1,Sibling!="!",(Mother=="!"&Father!="!")|(Father=="!"&Mother!="!")) %>% group_by(SYMBOL) %>% 
  arrange(IMPACT,SYMBOL)

final_table2<-final_table %>% group_by(SYMBOL) %>% summarise(variants_on_gene=n()) %>% left_join(final_table) %>% 
  dplyr::select(SYMBOL,variants_on_gene,Proband,Sibling,Mother,Father,
                Consequence, BIOTYPE, IMPACT,Feature, Feature_type,
                chr, pos, ref, alt, QUAL, Alleles, HGVS.Simple, HGVS.Protien,Amino_acids,
                SYMBOL, Gene,  EXON, INTRON, cDNA_position,
                CDS_position, Protein_position,  Codons, Existing_variation, DISTANCE, STRAND,     
                FLAGS, VARIANT_CLASS, SYMBOL_SOURCE, HGNC_ID, CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL,  
                TSL, APPRIS, CCDS, ENSP, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, GENE_PHENO,
                SIFT, PolyPhen, DOMAINS, miRNA, AF,
                gnomAD_AF, MAX_AF, MAX_AF_POPS, CLIN_SIG,             
                SOMATIC, PHENO, PUBMED, HGMD_link,gnomad_links,db_best_guess,omim_link,pubmed_links,genecard_links)


subset<-subset %>%dplyr::select(SYMBOL,variants_on_gene,Proband,Sibling,Mother,Father,
                              Consequence, BIOTYPE, IMPACT,Feature, Feature_type,
                        chr, pos, ref, alt, QUAL, Alleles, HGVS.Simple, HGVS.Protien,
                               SYMBOL, Gene,Amino_acids,  EXON, INTRON, cDNA_position,
                              CDS_position, Protein_position,  Codons, Existing_variation, DISTANCE, STRAND,     
                              FLAGS, VARIANT_CLASS, SYMBOL_SOURCE, HGNC_ID, CANONICAL, MANE_SELECT, MANE_PLUS_CLINICAL,  
                              TSL, APPRIS, CCDS, ENSP, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, GENE_PHENO,
                              SIFT, PolyPhen, DOMAINS, miRNA, AF,
                              gnomAD_AF, MAX_AF, MAX_AF_POPS, CLIN_SIG,             
                              SOMATIC, PHENO, PUBMED,HGMD_link,gnomad_links,db_best_guess,omim_link,pubmed_links,genecard_links)

write.csv(final_table2,"family18712_7-14-2022.csv",row.names = F)
write.csv(subset,"subset_family18712_7-7-2022.csv",row.names = F)
final_table2 %>% filter(SYMBOL%in%c("C20orf96","CES1","DDX60L","RDX","SPTBN5")) %>% 
  write.csv("family18712_SELECT_GENES.csv",row.names = F)
