library(shiny)
library(tidyverse)
library(igvShiny)
library(GenomicAlignments)
library(VariantAnnotation)
#library(igvR)
#-------------------------------------------------------------------------------
if(!dir.exists("tracks"))
    dir.create("tracks")
addResourcePath("tracks", "tracks")
#-------------------------------------------------------------------------------
f <- system.file(package="igvShiny")
stopifnot(file.exists(f))
printf <- function(...) print(noquote(sprintf(...)))

 
# Define UI for pathogenic variant table display
ui<-fluidPage(
    titlePanel("Variants of Interest"),
    
    # Create a new row for the table.
    tabsetPanel(
        tabPanel("Main App",
                 mainPanel(
                     fluidRow(
                         textOutput("Description"),
                     actionButton("addLocalVCFTrackButton", "Add VCF Trio"),
                     actionButton("addBamLocalFileButton", "Add Bam for Proband"),
                     actionButton("addBamLocalFileButton2", "Add Bam for Father"),
                     actionButton("addBamLocalFileButton3", "Add Bam for Mother"),
                     actionButton("removeUserTracksButton", "Remove User Tracks"),
                     igvShinyOutput('igvShiny_1'),
                     width=40),
                 fluidRow(checkboxInput("variable1", "Display Amino Acid and Transcript Information (ref and alt amino acid, cDNA and codon)",FALSE),
                          checkboxInput("variable2", "Display Gene info and Trio status",FALSE),
                          checkboxInput("variable3", "Display Genetic Algorthms (phasCons, pLI, oe, phyloP, cadd)",FALSE),
                          checkboxInput("variable4", "Display Allele Frequency Information",FALSE),
                          checkboxInput("variable5", "Display Refseq, Clinvar, and other links",FALSE)),
                 fluidRow(DT::dataTableOutput("variant_table")))),
        
        tabPanel("Data Dictionary",
                 fluidRow(tableOutput("dict"))),
        tabPanel("About",
                 mainPanel(
                 fluidRow(textOutput("more_about")),
                 fluidRow(imageOutput("About_img")))
        )
        )
)
# Define server logic required to create DT
#-------------------------------------------------------------------------------
server <- function(input, output, session) {
   output$Description<-renderText("This web app is designed for the exploration of possibly pathogenic variants. The variants are identified by a filtering strategy described on the ‘about’ tab. The application consists of an IGV instance and a table of variants. IGV stands for Integrative Genomics Viewer (igv.org). Currently it can be used to view VCF or BED files (it may be best to load one at a time due to memory constraints). The table is located below the viewer. Toggle the checks to see different information about each variant. Field descriptions are available on the ‘dictionary’ tab. The table is searchable with a global search field and each column can be sorted. Additionally, clicking on a table entry will cause the IGV instance to focus on the selected variant. ")
     observeEvent(input$addLocalVCFTrackButton, { 
        vcf_1<-readVcf("allfam.vcf","hg19")
        loadVcfTrack(session,id="igvShiny_1",trackName = "Trio",vcfData=vcf_1)
    })
    
    observeEvent(input$addBamLocalFileButton, {
        bam_1<-readGAlignments("2002145.bam")
        loadBamTrackFromLocalData(session,id="igvShiny_1",trackName = "Proband",data = bam_1)
    })
    
    observeEvent(input$addBamLocalFileButton2, {
        bam_2<-readGAlignments("2002161.bam")
        loadBamTrackFromLocalData(session,id="igvShiny_1",trackName = "Father",data = bam_2)
    })
    
    observeEvent(input$addBamLocalFileButton3, {
        bam_3<-readGAlignments("2002155.bam")
        loadBamTrackFromLocalData(session,id="igvShiny_1",trackName = "Mother",data = bam_3)
    })
    
    observeEvent(input$variant_table_rows_selected, {
        bam_1<-paste(variants[input$variant_table_rows_selected,1],variants[input$variant_table_rows_selected,2],sep = ":")
        showGenomicRegion(session,id="igvShiny_1",bam_1)
    })
    
    observeEvent(input$removeUserTracksButton, {
        printf("---- removeUserTracks")
        removeUserAddedTracks(session, id="igvShiny_1")
    })
    
    variants<-read_csv("pathogenic_variants_report.csv")
    dictionary<-read_csv("Data_Dictionary.csv",na = character())
    ID<-names(variants)[1:7]
    Amino_acid<-names(variants)[c(8:11,26:28)]
    trio<-names(variants)[12:18]
    algos<-names(variants)[19:23]
    gnom_afs<-names(variants)[24:25]
    the_rest<-names(variants)[29:46]
    var_table<-reactive({x<-ID
        if(input$variable1==T){x<-append(x,Amino_acid)}
        if(input$variable2==T){x<-append(x,trio)}
        if(input$variable3==T){x<-append(x,algos)}
        if(input$variable4==T){x<-append(x,gnom_afs)}
        if(input$variable5==T){x<-append(x,the_rest)}
        as.data.frame(variants[,x])
        })
    
    output$variant_table<-DT::renderDataTable(DT::datatable(data = var_table(),rownames = F,selection = "single",options = list(pageLength=5)))
    output$igvShiny_1<-renderIgvShiny(igvShiny(list(genomeName="hg19",initialLocus="chr1:1-50000")))
    output$dict<-renderTable(dictionary)
    output$more_about<-renderText("\n
    Whole exome sequencing data was obtained from a patient and their parents in a clinical setting. This data was made available in the form of compressed reference alignment map (CRAM) files. Using a bioinformatics package, Samtools, the data was converted into a variant call format (VCF) files and using the bioinformatics tool bcftools, merged into a single file. A VCF file only stores variations from a reference genome rather than an entire sequence. 
	The VCF file was than annotated. Annotation is the process of adding known biological information about the genome to the target data set. The VCF file was annotated using Bystro and the Variant Effect Predictor (VEP) tool from the Ensembl Genomes Project. The annotation processes pulls data from the NCBI Reference Sequence (refSeq) project, the Single Nucleotide Polymorphism Database (dbSNP),  the Clinvar database, and the genome Aggregation Database (gnomAD). 
	The refSeq database gives information about the structural effect of variants and the location in the genome. The dbSNP database gives some overlapping information with refSeq but also gives important single nucleotide polymorphisms including unique identifiers. Clinvar contains important information on disease causing or pathogenic variants. GnomAD contains important information on allele frequency in different populations.
	After annotation, the data was filtered to highlight variants that are more likely to be pathogenic. The filtering scheme was based on existing literature and was separated into primary and secondary filtering. Variants were filtered based on allele frequency, VCF quality scores, family segregation, cadd scores, phastCons score, phyloP scores, allele function, allele site type, and the variants presence in the clinvar database.
	\n
	\n
	The first filter was removing variants below a certain quality threshold. The VCF quality scores were base on a PHRED scaling. PHRED quality is equal to -10 log error probability so a higher quality implies a lower error probability. A quality score of 10 corresponds to 90% accuracy and  20  to 99%. Variants with a quality score below 20 were excluded. Family segregation means Variants were excluded if they were not present in the target individual, only present in their parents. Additionally, any variant with an clinvar entry with a clinical significance rating other than “benign” was automatically included in the final outcome. Variants could exist on multiple transcripts and therefore could have multiple values for certain fields. Variants were excluded if they did not have a listed function in refSeq and were only located on intergenic sites. Novel variants or variants with allele frequency of less than 1% were included in the final output table. 
	For secondary filtering, three different genome wide variant scores were used; cadd, phastCons, and phyloP. Cadd is the Combined Annotation Dependent Depletion score for determining the deleteriousness of simple insertions, deletions, or single nucleotide variants. Deleterious mutations are mutations that increase the risk of disease. Cadd scores are scaled similarly to PHRED scores. A cadd score of 20 implies a variant is in the 1% most deleterious variants and 30 would be among the 0.1%. A cutoff of 15 was used inline with recommendations from cadd’s creators, but pathogenic cutoffs are understood to be fairly arbitrary. Both phastCons and phyloP are algorithms that give conservation scores (Ramani 2019). Conservation scores try to measure how quickly or slowly a region of the genome evolves, with the assumption that slowly evolving regions are evolutionary important and variants in these regions could be harmful. PhastCons accounts for the effects of neighboring bases while phyloP does not. PhyloP is scaled 0 to 1 with values nearer to 1 evolving slowly. Phylop is scaled from -20 to 30 with positive values evolving slower than expected and negative values  evolving faster. Variants were included if they had phyloP scores above 3 or PhastCons scores greater than or equal to 0.9. A variant was filtered out only if was below all three score thresholds. 
	After filtering, some additional cleaning and merges took place. This could be considered an additional annotation step. Gnomad Loss of Function (lof) tables are designed to give additional information about frameshift, splice donor, splice acceptor, and stop-gain variants. Specifically these tables aim to give more information about whether variants cause loss of function for the proteins they code for. These tables added probability of loss of function intolerance (pLI) scores and observed expectation ratios (oe). PLI scores attempt to identify genes that are cannot tolerate truncating mutations. A gene with that changes phenotypes after a single loss of function mutation is known as a haploinsufficient gene, and pLI can be interpreted as the probability of a gene being haploinsufficient.  Oe ratios show the difference in the amount of observed LoF variants to the expected amount of LoF varaints if the rate of LoF variants was solely governed by chance.
	Information from the Human Gene Nomenclature Committee was added to get gene descriptions and allow links to be created to outside sources. Links were created to the Gene Cards website (https://www.genecards.org), pubmed (https://pubmed.ncbi.nlm.nih.gov), the Online Mendelian Inheritance in Man data base (OMIM)(https://www.omim.org), the Human Gene Mutation Database(HGMD)(www.hgmd.cf.ac.uk),the dbsnp database (www.ncbmi.nlm.nig.gov/snp), and gnomAD (https://gnomad.broadinstitute.org).
")
    output$About_img<-renderImage({
        filename<-"about.jpg"
        list(src=filename,
             alt="Data pipeline")},
        deleteFile=FALSE
    )
}


# Run the application 
shinyApp(ui = ui, server = server)
