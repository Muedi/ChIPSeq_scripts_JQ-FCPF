
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
# for (i in 1){

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
      dbobj <- dba.contrast(dbobj, minMembers=2)
      contrasts <- dba.show(dbobj, bContrasts = T)
      # run analysis
      dbobj <- dba.analyze(dbobj)
      dba.show(dbobj, bContrasts=TRUE)

      diffBound <- dba.report(dbobj)
      diffBound <- addGeneIDs(diffBound, "org.Hs.eg.db", "symbol")

      diffBound <- data.frame(diffBound)
      write.csv(paste0(out, "differentially_bound.csv"))
}
