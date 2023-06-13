######################################################################################################
# Script: ChipSeq Annotation
# Author: Maximilian Sprang, Muedi
# Date: 23.02.2023
# Description: This script annotatwes the Peakdata and produces exploratory plots, such as TSS profiles.
# Subsequently peaks are checked for enrichment in Biological Processes and Molecular function with the GO databse
######################################################################################################

library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
library(clusterProfiler)

# Specify the transcript database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Specify the folder containing the peak files
data.folder <- "peaks"

# Specify the genes of special interest
gene_list <- c("CMYC", "MYC", "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17")

# Get the peak files
peakFiles <- list.files(file.path(data.folder))
peakFiles <- file.path(data.folder, peakFiles)
#################################################################################################
# Annotation
#################################################################################################

# Perform annotation for each peak file
ann <- lapply(peakFiles, annotatePeak, TxDb = txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", tssRegion = c(-3000,3000)) 
names(ann) <- gsub(".bed$", "", basename(peakFiles))

# Subset annotation with genes of interest
ann_sub_list <- vector("list", length(names(ann)))
names(ann_sub_list) <- names(ann)
for (i in 1:length(names(ann))){
  ann_sub_list[[i]] <- ann[[i]]@anno[ann[[i]]@anno$SYMBOL %in% gene_list]
}
ann_sub_list <- Filter(function(x) length(x) > 0, ann_sub_list)
ann_sub <- lapply(ann_sub_list, annotatePeak, TxDb = txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", tssRegion = c(-3000,3000)) 

# Generate exploratory plots
# Plot feature distribution
plt <- plotAnnoBar(ann)
ggsave("output/annotation_bar.pdf", plot=plt, width=4, height=4)

# Plot distribution of TF-binding loci relative to TSS
plt <- plotDistToTSS(ann)
ggsave("output/dist_tss_bar.pdf", plot=plt, width=4, height=4)

# Calculate read count frequency with respect to the TSS
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(peakFiles, getTagMatrix, windows = prom)
names(tagMatrixList) <- gsub(".bed$", "", basename(peakFiles))
names(tagMatrixList)

# Generate average profile plot for all tag matrices
plot <- plotAvgProf(tagMatrixList, xlim = c(-3000, 3000))
ggsave("output/avgProf_TSS.pdf", plot=plot, width=10, height=5, device="pdf")

# Generate average profile plot for specific tag matrices
plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "3a_peaks", "4b_peaks")], xlim = c(-3000, 3000))
ggsave("output/avgProf_TSS_1.pdf", plot=plot, width=10, height=5, device="pdf")

plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "2a_peaks", "5b_peaks")], xlim = c(-3000, 3000))
ggsave("output/avgProf_TSS_2.pdf", plot=plot, width=10, height=5, device="pdf")

# Plot tag heatmap plots for specific tag matrices
pdf("output/tagMatrix_HM_1.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("1a_peaks", "1b_peaks")], xlim = c(-3000, 3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_2.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("2a_peaks", "2b_peaks")], xlim = c(-3000, 3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_3.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("3a_peaks", "3b_peaks")], xlim = c(-3000, 3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_4.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("4a_peaks", "4b_peaks")], xlim = c(-3000, 3000), color = "darkorange")
dev.off()

pdf("output/tagMatrix_HM_5.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("5a_peaks", "5b_peaks")], xlim = c(-3000, 3000), color = "darkorange")
dev.off()


#################################################################################################
# Enrichment
#################################################################################################

# Initialize lists to store enrichment results
BP_list <- list()
MF_list <- list()
peaknames <- names(ann)

# Perform enrichment analysis for each peak
for (i in 1:length(peaknames)) {
   # Extract gene annotations for the current peak
   anno <- data.frame(ann[[peaknames[i]]]@anno)
   entrez <- anno$geneId

   # Perform Gene Ontology (GO) enrichment analysis for Biological Process (BP) ontology
   ego <- enrichGO(gene = entrez,
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

   # Store the enrichment results in BP_list
   BP_list[peaknames[i]] <- ego
   cluster_summary <- data.frame(ego)
   write.csv(cluster_summary, paste0("output/clusterProfiler_BP_", peaknames[i], ".csv"))

   plt <- dotplot(ego, x = "count", showCategory = 20, font.size = 10)
   ggsave(paste0("output/clusterProfiler_BP_", peaknames[i], ".pdf"), plot = plt, dpi = 300, width = 6, height = 7)

   # Perform GO enrichment analysis for Molecular Function (MF) ontology
   ego <- enrichGO(gene = entrez,
                   keyType = "ENTREZID",
                   OrgDb = org.Hs.eg.db,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable = TRUE)

   # Store the enrichment results in MF_list
   MF_list[peaknames[i]] <- ego
   cluster_summary <- data.frame(ego)
   write.csv(cluster_summary, paste0("output/clusterProfiler_MF_", peaknames[i], ".csv"))

   plt <- dotplot(ego, x = "count", showCategory = 20, font.size = 10)
   ggsave(paste0("output/clusterProfiler_MF_", peaknames[i], ".pdf"), plot = plt, dpi = 300, width = 6, height = 7)
}

# Generate panel Dotplots for BP_list and MF_list
library(gridExtra)

# BP_list comparisons
pltlist <- BP_list[c("1b_peaks", "3a_peaks", "4b_peaks")]
grobs <- lapply(pltlist,
                function(.x) dotplot(.x, x = "count", showCategory = 20, font.size = 5) +
                  theme(legend.text = element_text(size = 5),
                        legend.title = element_text(size = 6),
                        legend.key.size = unit(0.3, 'cm')))

grid <- grid.arrange(grobs = grobs, ncol = 3)
path <- file.path("output/grid.GO.BP.1.pdf")
ggsave(path, grid, width = 10, height = 5, dpi = 300, device = 'pdf')

pltlist <- BP_list[c("1b_peaks", "2a_peaks", "5b_peaks")]
grobs <- lapply(pltlist,
                function(.x) dotplot(.x, x = "count", showCategory = 20, font.size = 5) +
                  theme(legend.text = element_text(size = 5),
                        legend.title = element_text(size = 6),
                        legend.key.size = unit(0.3, 'cm')))

grid <- grid.arrange(grobs = grobs, ncol = 3)
path <- file.path("output/grid.GO.BP.2.pdf")
ggsave(path, grid, width = 10, height = 5, dpi = 300, device = 'pdf')

# MF_list comparisons
pltlist <- MF_list[c("1b_peaks", "3a_peaks", "4b_peaks")]
grobs <- lapply(pltlist,
                function(.x) dotplot(.x, x = "count", showCategory = 20, font.size = 5) +
                  theme(legend.text = element_text(size = 5),
                        legend.title = element_text(size = 6),
                        legend.key.size = unit(0.3, 'cm')))

grid <- grid.arrange(grobs = grobs, ncol = 3)
path <- file.path("output/grid.GO.MF.1.pdf")
ggsave(path, grid, width = 10, height = 5, dpi = 300, device = 'pdf')

pltlist <- MF_list[c("1b_peaks", "2a_peaks", "5b_peaks")]
grobs <- lapply(pltlist,
                function(.x) dotplot(.x, x = "count", showCategory = 20, font.size = 5) +
                  theme(legend.text = element_text(size = 5),
                        legend.title = element_text(size = 6),
                        legend.key.size = unit(0.3, 'cm')))

grid <- grid.arrange(grobs = grobs, ncol = 3)
path <- file.path("output/grid.GO.MF.2.pdf")
ggsave(path, grid, width = 10, height = 5, dpi = 300, device = 'pdf')


#################################################################################################
# subset
# Generate the same plots for the subset of peaks for genes of interest
#################################################################################################

peaknames <- names(ann_sub)
# Feature distribution
plt <- plotAnnoBar(ann_sub)
ggsave("output/genes_of_interest/annotation_bar.pdf",plot=plt, width=4, height=4)
# Distribution of TF-binding loci relative to TSS
plt <- plotDistToTSS(ann_sub)
ggsave("output/genes_of_interest/dist_tss_bar.pdf", plot=plt, width=4, height=4)


# Calculate read count frequency with respect to the TSS
prom <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)
tagMatrixList <- lapply(ann_sub_list, getTagMatrix, windows = prom)
 tagMatrixList <- tagMatrixList[2:8] # drop length 0

# Generate average profile plot for all tag matrices
plot <- plotAvgProf(tagMatrixList, xlim = c(-3000,3000))
ggsave("output/genes_of_interest/avgProf_TSS.pdf", plot=plot, width=10, height=5, device="pdf")

# Generate average profile plot for specific tag matrices
plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "3a_peaks", "4b_peaks")], xlim = c(-3000,3000))
ggsave("output/genes_of_interest/avgProf_TSS_1.pdf", plot=plot, width=10, height=5, device="pdf")

plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "2a_peaks", "5b_peaks")], xlim = c(-3000,3000))
ggsave("output/genes_of_interest/avgProf_TSS_2.pdf", plot=plot, width=10, height=5, device="pdf")


# Plot tag heatmap plots for specific tag matrices
pdf("output/genes_of_interest/tagMatrix_HM_1.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("1a_peaks", "1b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/genes_of_interest/tagMatrix_HM_2.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("2a_peaks", "2b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/genes_of_interest/tagMatrix_HM_3.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("3a_peaks", "3b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/genes_of_interest/tagMatrix_HM_4.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("4a_peaks", "4b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()

pdf("output/genes_of_interest/tagMatrix_HM_5.pdf", width=10, height=5)
tagHeatmap(tagMatrixList[c("5a_peaks", "5b_peaks")], xlim = c(-3000,3000), color = "darkorange")
dev.off()


#################################################################################################
# Enrichment
#################################################################################################
BP_list <- list()
MF_list <- list()

for (i in 1:length(peaknames)) {
   anno <- data.frame(ann_sub[[peaknames[i]]]@anno)
   entrez <- anno$geneId

   ego <- enrichGO(gene = entrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

    BP_list[peaknames[i]] <- ego
    cluster_summary <- data.frame(ego)   
    write.csv(cluster_summary, paste0("output/genes_of_interest/clusterProfiler_BP_", peaknames[i], ".csv"))        

    plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
    ggsave(paste0("output/genes_of_interest/clusterProfiler_BP_", peaknames[i], ".pdf"), plot=plt, dpi=300, width=6, height=7)


   ego <- enrichGO(gene = entrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

    MF_list[peaknames[i]] <- ego
    cluster_summary <- data.frame(ego)   
    write.csv(cluster_summary, paste0("output/genes_of_interest/clusterProfiler_MF_", peaknames[i], ".csv"))        

    plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
    ggsave(paste0("output/genes_of_interest/clusterProfiler_MF_", peaknames[i], ".pdf"), plot=plt, dpi=300, width=6, height=7)

}

library(gridExtra)
# comparisons
pltlist <- BP_list[c("1b_peaks", "3a_peaks", "4b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/genes_of_interest/grid.GO.BP.1.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- BP_list[c("1b_peaks", "2a_peaks", "5b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/genes_of_interest/grid.GO.BP.2.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- MF_list[c("1b_peaks", "3a_peaks", "4b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/genes_of_interest/grid.GO.MF.1.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- MF_list[c("1b_peaks", "2a_peaks", "5b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/genes_of_interest/grid.GO.MF.2.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')
