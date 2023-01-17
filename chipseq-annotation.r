#### chipseq annotation
library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Mm.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

data.folder <- "peaks/"

# get files
files <- 

chip_peaks_ranges <- toGRanges(file.path(data.folder, ))
chip_peaks_seeker <- annotatePeak(chip_peaks_ranges, 
                                  TxDb = txdb, 
                                  genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic"), 
                                  annoDb = "org.Mm.eg.db", 
                                  tssRegion = c(-3000,3000))