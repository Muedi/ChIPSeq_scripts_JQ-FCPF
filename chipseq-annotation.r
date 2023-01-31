#### chipseq annotation
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

data.folder <- "peaks"

# get files
peakFiles <- list.files(file.path(data.folder))
peakFiles <- file.path(data.folder, peakFiles)
#################################################################################################
# annotation
#################################################################################################
ann <- lapply(peakFiles, annotatePeak, TxDb = txdb, genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), annoDb = "org.Hs.eg.db", tssRegion = c(-3000,3000)) 
names(ann) <- gsub(".bed$", "", basename(peakFiles))
peaknames <- names(ann)
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

plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "3a_peaks", "4b_peaks")], xlim = c(-3000,3000))
ggsave("output/avgProf_TSS_1.pdf", plot=plot, width=10, height=5, device="pdf")

plot <- plotAvgProf(tagMatrixList[c("1b_peaks", "2a_peaks", "5b_peaks")], xlim = c(-3000,3000))
ggsave("output/avgProf_TSS_2.pdf", plot=plot, width=10, height=5, device="pdf")



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


#################################################################################################
# enrichment
#################################################################################################
library(clusterProfiler)
library(org.Hs.eg.db)

BP_list <- list()
MF_list <- list()

for (i in 1:length(peaknames)) {
# for (i in 1) {
   anno <- data.frame(ann[[peaknames[i]]]@anno)
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
    write.csv(cluster_summary, paste0("output/clusterProfiler_BP_", peaknames[i], ".csv"))        

    plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
    ggsave(paste0("output/clusterProfiler_BP_", peaknames[i], ".pdf"), plot=plt, dpi=300, width=6, height=7)


   ego <- enrichGO(gene = entrez, 
                    keyType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db, 
                    ont = "MF", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE)

    MF_list[peaknames[i]] <- ego
    cluster_summary <- data.frame(ego)   
    write.csv(cluster_summary, paste0("output/clusterProfiler_MF_", peaknames[i], ".csv"))        

    plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
    ggsave(paste0("output/clusterProfiler_MF_", peaknames[i], ".pdf"), plot=plt, dpi=300, width=6, height=7)

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
path <- file.path("output/grid.GO.BP.1.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- BP_list[c("1b_peaks", "2a_peaks", "5b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/grid.GO.BP.2.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- MF_list[c("1b_peaks", "3a_peaks", "4b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/grid.GO.MF.1.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

pltlist <- MF_list[c("1b_peaks", "2a_peaks", "5b_peaks")]

grobs <- lapply(pltlist,
        function(.x) dotplot(.x, x="count", showCategory=20, font.size=5) + 
        theme(legend.text=element_text(size=5),
              legend.title=element_text(size=6),
              legend.key.size = unit(0.3, 'cm')) )

grid <- grid.arrange(grobs=grobs, ncol=3) 
path <- file.path("output/grid.GO.MF.2.pdf")
ggsave(path, grid, width=10, height=5,  dpi = 300, device='pdf')

#################################################################################################
# Differential peaks
#################################################################################################
library(DiffBind)
library(stringr)


out.dir <- "output/diff"
dir.create(out.dir)

# samples <- read.csv("samplesheets/comp/samplesheet_all.csv")
# samples_no1a <- read.csv("samplesheets/comp/samplesheet_-1a.csv")
# get sample sheets 
sheets <- list.files(file.path("samplesheets"))
sheets <- file.path("samplesheets", sheets)
sheets <- sheets[!file.info(sheets)$isdir]
names(sheets) <- str_extract_all(sheets, "[0-9]v[0-9]")


for (i in 1:length(sheets)){

      # make output
      out <- paste0(out.dir, names(sheets[i]))
      samples <- read.csv(sheets[[i]])
      # init obj
      dbobj <- dba(sampleSheet=samples)
      dbobj$config$yieldSize <- 500000
      dbobj$config$cores <- 3
      # corr 
      pdf(paste0(out, "correlatHM.pdf"), width=7, height=7)
      plot(dbobj)
      dev.off()
      # count
      dbobj <- dba.count(dbobj)
      # plot from readcounts rather than peaks:
      pdf(paste0(out, "correlatHM_readCounts.pdf"), width=7, height=7)
      plot(dbobj)
      dev.off()
      # save information
      info <- dba.show(dbobj)
      # normalize
      dbobj <- dba.normalize(dbobj)
      # contrasts
      dbobj <- dba.contrst(dbobj, minMembers=1)
      contrasts <- dba.show(diff_WThox_KOv2, bContrasts = T)
      # run analysis
      dbobj <- dba.analyze(dbobj)
      dba.show(dbobj, bContrasts=TRUE)

      diffBound <- dba.report(dbobj)
      diffBound <- addGeneIDs(diffBound, "org.Hs.eg.db", "symbol")

      diffBound <- data.frame(diffBound)
      write.csv(paste0(out, "differentially_bound.csv"))
}
