######################################################################################################
# Script: ChipSeq Combi plts HM and AvgProf
# Author: Maximilian Sprang, Muedi
# Date: 23.02.2023
# Description: This script produces plots that combine the TSS Heatmaps and the Average Peak Profile around the TSS.
######################################################################################################

# Load required libraries
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
# library(clusterProfiler)
library(gridExtra)

# Set the reference annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Set the folder containing the peak files
data.folder <- "peaks"

# Get the list of peak files in the data folder
peakFiles <- list.files(file.path(data.folder))
peakFiles <- file.path(data.folder, peakFiles)

# Read count frequency with respect to the TSS
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(peakFiles, getTagMatrix, windows = prom)
names(tagMatrixList) <- gsub(".bed$", "", basename(peakFiles))
peaknames <- names(tagMatrixList)

# Loop through each peak file
for (i in 1:length(peaknames)) {
    
    inp <- tagMatrixList[peaknames[i]]
    names(inp) <- NULL

    # Generate a heatmap of the peak intensities
    hm <- as.grob( 
        function()   tagHeatmap(inp, xlim = c(-3000,3000), color = "darkorange")
    )
    # Generate an average profile plot of the peak intensities
    avg <- plotAvgProf(tagMatrixList[peaknames[i]], xlim = c(-3000,3000), color = "darksteelblue") + 
        theme(legend.position="none")
    
    # Arrange the heatmap and average profile plot side by side
    out <- grid.arrange(grobs = list(avg, hm), ncol=2, top=peaknames[i])

    # Save the plot as PNG and PDF files
    ggsave(paste0("output/avg_hm_plots/", peaknames[i], ".png"), plot=out, width=10, height=4, dpi=300)
    ggsave(paste0("output/avg_hm_plots/", peaknames[i], ".pdf"), plot=out, width=10, height=4, dpi=300)
}
