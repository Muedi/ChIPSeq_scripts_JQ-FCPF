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

# genes of special interest
gene_list <- c( "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17")
gene_entrz <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = gene_list,
                                    columns = "ENTREZID",
                                    keytype="SYMBOL")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak.folder <- "peaks"
bam.folder <- "bams"

gtrack <- GenomeAxisTrack()

txdb_gene_interest <- AnnotationDbi::select(txdb, 
                                            keys = gene_entrz$ENTREZID,
                                            columns=columns(txdb),
                                            keytype="GENEID")
txdb_gene_interest <- txdb_gene_interest %>% as_tibble()

for (i in 1:length(gene_list)) {
    symbol <- gene_list[i]
    id <- gene_entrz[gene_entrz$SYMBOL == symbol,"ENTREZID"]

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
        unique()
    end <- txdb_gene_interest %>%
        dplyr::filter(GENEID == id) %>%
        pull("EXONEND") %>%
        unique()
    txTr <- GeneRegionTrack(txdb,
                            chromosome = chr, 
                            start = min(start),
                            end = max(end))
    plotTracks(c(gtrack, txTr))
}



plotTracks(c(gtrack, txTr))
