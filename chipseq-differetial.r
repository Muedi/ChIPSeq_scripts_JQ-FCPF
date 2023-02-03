
#################################################################################################
# Differential peaks
#################################################################################################
library(DiffBind)
library(stringr)
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(ggplot2)


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

out.dir <- "output/diff/"
dir.create(out.dir)

# samples <- read.csv("samplesheets/comp/samplesheet_all.csv")
# samples_no1a <- read.csv("samplesheets/comp/samplesheet_-1a.csv")
# get sample sheets 
sheets <- list.files(file.path("samplesheets"))
sheets <- file.path("samplesheets", sheets)
sheets <- sheets[!file.info(sheets)$isdir]
names(sheets) <- str_extract_all(sheets, "[0-9]v[0-9]")


for (i in 2:length(sheets)){
# for (i in 1){

      # make output
      out <- paste0(out.dir, names(sheets[i]))
      samples <- read.csv(sheets[[i]])
      # init obj
      dbobj <- dba(sampleSheet=samples)
      dbobj$config$yieldSize <- 400000
      dbobj$config$cores <- 3
      # corr 
      pdf(paste0(out, "_correlatHM.pdf"), width=7, height=7)
      plot(dbobj)
      dev.off()
      # count
      dbobj <- dba.count(dbobj)
      # plot from readcounts rather than peaks:
      pdf(paste0(out, "_correlatHM_readCounts.pdf"), width=7, height=7)
      plot(dbobj)
      dev.off()
      # save information
      info <- dba.show(dbobj)
      # normalize
      dbobj <- dba.normalize(dbobj)
      # contrasts
      dbobj <- dba.contrast(dbobj, minMembers=2)
      contrasts <- dba.show(dbobj, bContrasts = T)
      # run analysis
      dbobj <- dba.analyze(dbobj)
      dba.show(dbobj, bContrasts=TRUE)

      diffBound <- dba.report(dbobj, contrast=1, th=0.1, bUsePval=T)
      diffBound <- annotatePeak(diffBound, TxDb = txdb,
            genomicAnnotationPriority = c("Promoter",
                                          "5UTR",
                                          "3UTR",
                                          "Exon",
                                          "Intron",
                                          "Downstream",
                                          "Intergenic"),
            annoDb = "org.Hs.eg.db",
            tssRegion = c(-3000,3000)) 


      diffBound <- data.frame(diffBound)
      write.csv(diffBound, paste0(out, "_differentially_bound.csv"))

      entrez <- diffBound$geneId

      ego <- enrichGO(gene = entrez, 
                        keyType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db, 
                        ont = "BP", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
      if (dim(ego)[1] > 0) {

      # BP_list[peaknames[i]] <- ego
      cluster_summary <- data.frame(ego)   
      write.csv(cluster_summary, paste0(out, "clusterProfiler_BP_", ".csv"))        

      plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
      ggsave(paste0(out, "clusterProfiler_BP_", ".pdf"), plot=plt, dpi=300, width=6, height=7)
      } else {
            print(paste0(names(sheets[i]), ": no enrichment found"))
      }


      ego <- enrichGO(gene = entrez, 
                        keyType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db, 
                        ont = "MF", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable = TRUE)
      if (dim(ego)[1] > 0) {

      # MF_list[peaknames[i]] <- ego
      cluster_summary <- data.frame(ego)   
      write.csv(cluster_summary, paste0(out, "clusterProfiler_MF_", ".csv"))        

      plt <- dotplot(ego, x="count", showCategory=20, font.size=10, )
      ggsave(paste0(out, "clusterProfiler_MF_", ".pdf"), plot=plt, dpi=300, width=6, height=7)
      } else {
            print(paste0(names(sheets[i]), ": no enrichment found"))
      }

}
