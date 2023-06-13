######################################################################################################
# Script: ChipSeq Trackplots for multiple genes
# Author: Maximilian Sprang, Muedi
# Date: 23.02.2023
# Description: This script produces trackplots to observe the differences of Binding between the different sample treatments for multiple target genes
######################################################################################################

library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
library(gridExtra)
library(Gviz)
library(tidyverse)
library(stringr)
library(biomaRt)
library(trackViewer)

# genes of special interest
# Select genes of interest based on their symbols and retrieve their ENTREZID using org.Hs.eg.db
gene_entrz <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = c("MYC", "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17", "FOS"),
                                    columns = "ENTREZID",
                                    keytype = "SYMBOL")

gene_list <- list("MYC", "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17", "FOS")

# mart to get locations
# Initialize ensembl mart and dataset to retrieve gene locations
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# peak.folder <- "peaks"
bam.folder <- "bams"
# initialize Axis
axistrack <- GenomeAxisTrack()

# Retrieve gene annotations for the genes of interest using TxDb
txdb_gene_interest <- AnnotationDbi::select(txdb, 
                                            keys = gene_entrz$ENTREZID,
                                            columns = columns(txdb),
                                            keytype = "GENEID")
txdb_gene_interest <- txdb_gene_interest %>% as_tibble()

bwfiles <- list.files(bam.folder, pattern = ".bw$", full.names = TRUE)
names(bwfiles) <- str_extract(bwfiles, "[0-9][ab]")

bamfiles <- list.files(bam.folder, pattern = ".bam$", full.names = TRUE)
names(bamfiles) <- str_extract(bamfiles, "[0-9][ab]")
# coverage track from BigWig files
sample_cov <- lapply(bwfiles, AlignmentsTrack, isPaired = FALSE, ucscChromosomeNames = FALSE)

options(ucscChromosomeNames = FALSE)

viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(0.07, 0.04, 0.01, 0.03))
setTrackViewerStyleParam(viewerStyle, "xaxis", TRUE)

max_overall <- 0
max <- 215

# trackviewer
for (i in 1:length(gene_list)) {
    symbol <- gene_list[[i]]
    id <- get(symbol, org.Hs.egSYMBOL2EG)

    # get loc
    coords <- getBM(attributes = c('chromosome_name', 'start_position', 'end_position', 'strand'),
                    filters = c('hgnc_symbol'),
                    values = list(symbol),
                    mart = ensembl)
    chr <- paste0("chr", coords$chromosome_name)
    start <- coords$start_position - 10000
    end <- coords$end_position + 10000 
    strand <- coords$strand

    # range
    gr <- GRanges(chr, IRanges(start, end), strand = strand)

    # gene model and track
    trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                             org.Hs.eg.db,
                             gr = gr)
    gTrack <- geneTrack(id, TxDb.Hsapiens.UCSC.hg38.knownGene, symbol)[[1]]
    setTrackStyleParam(gTrack, "color", "darkblue")
    setTrackStyleParam(gTrack, "height", 0.001)
    setTrackStyleParam(gTrack, "ylabgp", list("cex" = 0.6, fontface = "bold"))

    # samples: (1b, 3a, 4b, 2b, and 5b)
    NT <- importScore(bwfiles["1b"], format = "BigWig", ranges = gr)
    setTrackStyleParam(NT, "color", "black")
    setTrackStyleParam(NT, "height", 0.1998)
    setTrackStyleParam(NT, "ylabgp", list("cex" = 0.6, fontface = "bold"))
    setTrackYaxisParam(NT, "label", TRUE)
    setTrackYaxisParam(NT, "gp", list("cex" = 0.6))

    JQ1_FCPF_low <- importScore(bwfiles["3a"], format = "BigWig", ranges = gr)
    setTrackStyleParam(JQ1_FCPF_low, "color", "darkorange")
    setTrackStyleParam(JQ1_FCPF_low, "height", 0.1998)
    setTrackStyleParam(JQ1_FCPF_low, "ylabgp", list("cex" = 0.6, fontface = "bold"))
    
    sgRNA_JQ1_FCPF_low <- importScore(bwfiles["4b"], format = "BigWig", ranges = gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "height", 0.1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "ylabgp", list("cex" = 0.6, fontface = "bold"))
    
    sgRNA_JQ1_FCPF_high <- importScore(bwfiles["2b"], format = "BigWig", ranges = gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "height", 0.1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "ylabgp", list("cex" = 0.6, fontface = "bold"))
    
    sgRNQ_JQ1_high <- importScore(bwfiles["5a"], format = "BigWig", ranges = gr)
    setTrackStyleParam(sgRNQ_JQ1_high, "color", "#9400D3")
    setTrackStyleParam(sgRNQ_JQ1_high, "height", 0.1998)
    setTrackStyleParam(sgRNQ_JQ1_high, "ylabgp", list("cex" = 0.6, fontface = "bold"))

    tl <- trackList(gTrack,
                    sgRNQ_JQ1_high,
                    sgRNA_JQ1_FCPF_high,
                    sgRNA_JQ1_FCPF_low,
                    JQ1_FCPF_low,
                    NT
    )
    
    names(tl) <- c(symbol,
                   "sgRNQ-JQ1 high",
                   "sgRNA-JQ1_FCPF high",
                   "sgRNA-JQ1_FCPF low",
                   "JQ1-FCPF low",
                   "NT"
    )
    # get max score
    # max <- 0
    # for(i in 1:length(tl)){
    #     #print(tl[[i]]$dat$score)
    #     new <- max(tl[[i]]$dat$score)
    #     if (new > max){
    #         max <- new
    #     }
    #     if (new > max_overall){
    #         max_overall <- new
    #     }
        
    # }

    # fixed y-axis
    for (j in 2:length(tl)) {
        setTrackStyleParam(tl[[j]], "ylim", c(0, max))
    }

    png(paste0("output/tracks/", symbol, ".png"), width = 3000, height = 2200, res = 300)
    viewTracks(tl, gr = gr, viewerStyle = viewerStyle)
    dev.off()

    svg(paste0("output/tracks/", symbol, ".svg"), width = 10, height = 7.3)
    viewTracks(tl, gr = gr, viewerStyle = viewerStyle)
    dev.off()
}

#    browseTracks(tl, gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)
