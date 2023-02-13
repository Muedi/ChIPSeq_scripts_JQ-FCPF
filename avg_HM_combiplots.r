#### chipseq annotation
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

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

data.folder <- "peaks"

# genes of special interest
gene_list <- c("CMYC", "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17")


# get files
peakFiles <- list.files(file.path(data.folder))
peakFiles <- file.path(data.folder, peakFiles)

# Read count frequency respect to the TSS
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(peakFiles, getTagMatrix, windows = prom)
names(tagMatrixList) <- gsub(".bed$", "", basename(peakFiles))
peaknames <- names(tagMatrixList)

for (i in 1:length(peaknames)) {
    
    inp <- tagMatrixList[peaknames[i]]
    names(inp) <- NULL

    hm <- as.grob( 
        function()   tagHeatmap(inp, xlim = c(-3000,3000), color = "darkorange")
    )
    avg <- plotAvgProf(tagMatrixList[peaknames[i]], xlim = c(-3000,3000), color = "darksteelblue") + 
        theme(legend.position="none")
    
    out <- grid.arrange(grobs = list(avg, hm), ncol=2, top=peaknames[i])
    ggsave(paste0("output/avg_hm_plots/", peaknames[i], ".png"), plot=out, width=10, height=4, dpi=300)
}