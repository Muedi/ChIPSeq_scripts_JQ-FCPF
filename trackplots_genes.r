library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library(gridExtra)
library(Gviz)
library(tidyverse)
library(stringr)
library(org.Hs.eg.db)

# genes of special interest
gene_list <- c( "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17", "FOS")
gene_entrz <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = gene_list,
                                    columns = "ENTREZID",
                                    keytype="SYMBOL")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak.folder <- "peaks"
bam.folder <- "bams"

axistrack <- GenomeAxisTrack()


txdb_gene_interest <- AnnotationDbi::select(txdb, 
                                            keys = gene_entrz$ENTREZID,
                                            columns=columns(txdb),
                                            keytype="GENEID")
txdb_gene_interest <- txdb_gene_interest %>% as_tibble()

bwfiles <- list.files("bams", pattern = ".bw$", full.names=T)
names(bwfiles) <- str_extract(bwfiles, "[0-9][ab]")

bamfiles <- list.files("bams", pattern = ".bam$", full.names=T)
names(bamfiles) <- str_extract(bamfiles, "[0-9][ab]")


sample_cov <- lapply(bwfiles, AlignmentsTrack, isPaired=F, ucscChromosomeNames=FALSE)

options(ucscChromosomeNames=FALSE)

# for (i in 1:length(gene_list)) {
#     symbol <- gene_list[i]
#     symbol <- "MYC"
#     id <- gene_entrz[gene_entrz$SYMBOL == symbol,"ENTREZID"]

#     # get adress
#     # chr should be singular after unqiue, while multiple starts/ends remain
#     # use the min start and the max end to cover the complete gene
#     chr <- txdb_gene_interest %>%
#         dplyr::filter(GENEID == id) %>%
#         pull("EXONCHROM") %>%
#         unique() 
#     start <- txdb_gene_interest %>%
#         dplyr::filter(GENEID == id) %>%
#         pull("EXONSTART") %>%
#         unique() %>% min() - 2000
#     end <- txdb_gene_interest %>%
#         dplyr::filter(GENEID == id) %>%
#         pull("EXONEND") %>%
#         unique() %>% max() + 2000
    
#     gtTrack <- GeneRegionTrack(txdb,
#                 chromosome=chr,
#                 start=min(start),
#                 end=max(end),
#                 collapseTranscripts="meta",
#                 transcriptAnnotation="name", 
#                 showId=TRUE)

#     knownGenes <- UcscTrack(genome = "hg38", chromosome = chr,
#                         track = "NCBI RefSeq", from = start, to = end,
#                         trackType = "GeneRegionTrack",
#                         rstarts = "exonStarts", rends = "exonEnds",
#                         gene = "name", symbol = "name",
#                         transcript = "name", strand = "strand",
#                         fill = "#8282d2", name = "UCSC Genes")


#     NT <- DataTrack(bwfiles["1b"], type="hist" , window = -1, windowSize = 250)
    
#     JQ1_high <- DataTrack(bwfiles["5a"], type="hist" , window = -1, windowSize = 250)

#     NT_algn <- AlignmentsTrack(bamfiles["1b"], 
#                                 chromosome=chr,
#                                 from = min(start),
#                                 to = max(end))
    
#     plt <- plotTracks(c(axistrack, NT, JQ1_high, knownGenes), from = min(start), to = max(end))


#     break
# }


# trackviewer
library(trackViewer)
for (i in 1:length(gene_list)) {
    symbol <- gene_list[i]
    symbol <- "MYC"
    id <- get(symbol, org.Hs.egSYMBOL2EG)

    # get adress
    # chr should be singular after unqiue, while multiple starts/ends remain
    # use the min start and the max end to cover the complete gene
    chr <- txdb_gene_interest %>%
        dplyr::filter(GENEID == id) %>%
        pull("EXONCHROM") %>%
        unique() 
    start <- txdb_gene_interest %>%
        dplyr::filter(GENEID == id) %>%
        pull("EXONSTART") %>%
        unique() %>% min() - 2000
    end <- txdb_gene_interest %>%
        dplyr::filter(GENEID == id) %>%
        pull("EXONEND") %>%
        unique() %>% max() + 2000
    
    gr = GRanges(chr, IRanges(start, end), strand="-")
    
    trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
    gTrack <- geneTrack(id,TxDb.Hsapiens.UCSC.hg38.knownGene, symbol)[[1]]

    NT <- importScore(bwfiles["1b"], format="BigWig", ranges=gr)
    
    JQ1_high <- importScore(bwfiles["5a"],  format="BigWig", ranges=gr)

    # NT_algn <- importBam(bamfiles["1b"])
    
    vp <- viewTracks(trackList(gTrack, NT, JQ1_high), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)
    svg("output/test.svg", width= 7, height = 5)
    print(vp)
    dev.off()
    break
}


