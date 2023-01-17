#### chipseq annotation
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Mm.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

data.folder <- "peaks"

# get files
peakFiles <- list.files(file.path(data.folder))
peakFiles <- file.path(data.folder, peakFiles)

ann <- lapply(peakFiles, annotatePeak, TxDb = txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", tssRegion = c(-3000,3000)) 
names(ann) <- gsub(".bed$", "", basename(peakFiles))
names(ann)
# Feature distribution
# upsetplot(ann[[1]])
plt <- plotAnnoBar(ann)
ggsave("output/annotation_bar.pdf",plot=plt, width=4, height=4)
# Distribution of TF-binding loci relative to TSS
plt <- plotDistToTSS(ann)
ggsave("output/dist_tss_bar.pdf", plot=plt, width=4, height=4)


# Read count frequency respect to the TSS
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(peakFiles, getTagMatrix, windows = prom)
names(tagMatrixList) <- gsub(".bed$", "", basename(peakFiles))
names(tagMatrixList)
 

plot <- plotAvgProf(tagMatrixList, xlim = c(-3000,3000))
ggsave("output/avgProf_TSS.pdf", plot=plot, width=10, height=5, device="pdf")


pdf("output/tagMatrix_HM_1.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("1a_peaks", "1b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_2.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("2a_peaks", "2b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_3.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("3a_peaks", "3b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_4.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("4a_peaks", "4b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_5.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("5a_peaks", "5b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()