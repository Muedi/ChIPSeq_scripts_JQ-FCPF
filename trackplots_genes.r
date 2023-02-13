library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)

library(Gviz)

# genes of special interest
gene_list <- c("CMYC", "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17")


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak.folder <- "peaks"
bam.folder <- "bams"

gtrack <- GenomeAxisTrack()

select(txdb, keys = gene_list, columns=)

txTr <- GeneRegionTrack(txdb, chromosome = "chr6", 
                        start = 35000000,  end = 40000000)

plotTracks(c(gtrack, txTr))
