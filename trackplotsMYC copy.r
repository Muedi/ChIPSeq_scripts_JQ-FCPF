library(ChIPseeker)
library(GenomicFeatures)
library(GenomicRanges)
library(org.Hs.eg.db)
library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(openxlsx)
library(ggplot2)
library(clusterProfiler)
library(gridExtra)
library(Gviz)
library(tidyverse)
library(stringr)
library(biomaRt)
library(trackViewer)

# genes of special interest
gene_entrz <- AnnotationDbi::select(org.Hs.eg.db, 
                                    keys = c( "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17", "FOS"),
                                    columns = "ENTREZID",
                                    keytype="SYMBOL")


gene_list <- list( "MYC",  "JUN", "YY1", "E2F3", "GSTO1", "NOP56", "NUDT17", "FOS")
# mart to  get locations
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak.folder <- "peaks"
bam.folder <- "bams"

axistrack <- GenomeAxisTrack()


txdb_gene_interest <- AnnotationDbi::select(txdb, 
                                            keys = gene_entrz$ENTREZID,
                                            columns=columns(txdb),
                                            keytype="GENEID")
txdb_gene_interest <- txdb_gene_interest %>% as_tibble()

bwfiles <- list.files("bams", pattern = ".bw$", full.names=T)
names(bwfiles) <- str_extract(bwfiles, "[0-9][ab]")

bamfiles <- list.files("bams", pattern = ".bam$", full.names=T)
names(bamfiles) <- str_extract(bamfiles, "[0-9][ab]")


sample_cov <- lapply(bwfiles, AlignmentsTrack, isPaired=F, ucscChromosomeNames=FALSE)

options(ucscChromosomeNames=FALSE)

viewerStyle <- trackViewerStyle()
setTrackViewerStyleParam(viewerStyle, "margin", c(.07, .04, .01, .03))
setTrackViewerStyleParam(viewerStyle, "xaxis", T)
# setTrackViewerStyleParam(viewerStyle, "autolas", T)


max_overall <- 0
max <- 215
# trackviewer
#for (i in 1:length(gene_list)) {
for (i in 1) {

    symbol <- gene_list[[i]]
    id <- get(symbol, org.Hs.egSYMBOL2EG)

    # get loc
    coords <- getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'strand'),
      filters=c('hgnc_symbol'),
      values=list(symbol),
      mart=ensembl)
    chr <- paste0("chr", coords$chromosome_name)
    start <- coords$start_position - 10000
    end <- coords$end_position + 10000 
    strand <- coords$strand


    # range
    gr = GRanges(chr, IRanges(start, end), strand=strand)
    # ideogram
    # ideo <- loadIdeogram("hg38", chrom=chr, ranges=IRanges(start, end))
    # gene model and track
    trs <- geneModelFromTxdb(TxDb.Hsapiens.UCSC.hg38.knownGene,
                         org.Hs.eg.db,
                         gr=gr)
    gTrack <- geneTrack(id,TxDb.Hsapiens.UCSC.hg38.knownGene, symbol)[[1]]
    setTrackStyleParam(gTrack, "color", "darkblue")
    setTrackStyleParam(gTrack, "height", .001)
    # setTrackStyleParam(gTrack, "ylabpos", "right")
    setTrackStyleParam(gTrack, "ylabgp", list("cex"=0.6, fontface="bold"))

    # samples: (1b, 3a, 4b, 2b and 5b
    NT <- importScore(bwfiles["1a"], format="BigWig", ranges=gr)
    setTrackStyleParam(NT, "color", "black")
    setTrackStyleParam(NT, "height", .1998)
    setTrackStyleParam(NT, "ylabgp", list("cex"=0.6, fontface="bold"))
    setTrackYaxisParam(NT, "label", T)
    setTrackYaxisParam(NT, "gp", list("cex"=0.6))

    NTb <- importScore(bwfiles["1b"], format="BigWig", ranges=gr)
    setTrackStyleParam(NTb, "color", "black")
    setTrackStyleParam(NTb, "height", .1998)
    setTrackStyleParam(NTb, "ylabgp", list("cex"=0.6, fontface="bold"))
    setTrackYaxisParam(NTb, "label", T)
    setTrackYaxisParam(NTb, "gp", list("cex"=0.6))

    JQ1_FCPF_low <- importScore(bwfiles["3a"], format="BigWig", ranges=gr)
    setTrackStyleParam(JQ1_FCPF_low, "color", "darkorange")
    setTrackStyleParam(JQ1_FCPF_low, "height", .1998)
    setTrackStyleParam(JQ1_FCPF_low, "ylabgp", list("cex"=0.6, fontface="bold"))
    
    JQ1_FCPF_low_b <- importScore(bwfiles["3b"], format="BigWig", ranges=gr)
    setTrackStyleParam(JQ1_FCPF_low_b, "color", "darkorange")
    setTrackStyleParam(JQ1_FCPF_low_b, "height", .1998)
    setTrackStyleParam(JQ1_FCPF_low_b, "ylabgp", list("cex"=0.6, fontface="bold"))
    
    sgRNA_JQ1_FCPF_low <- importScore(bwfiles["4a"], format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low, "ylabgp", list("cex"=0.6, fontface="bold"))

    sgRNA_JQ1_FCPF_low_b <- importScore(bwfiles["4b"], format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low_b, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_low_b, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_low_b, "ylabgp", list("cex"=0.6, fontface="bold"))
    
    sgRNA_JQ1_FCPF_high <- importScore(bwfiles["2a"], format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high, "ylabgp", list("cex"=0.6, fontface="bold"))
        
    sgRNA_JQ1_FCPF_high_b <- importScore(bwfiles["2b"], format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high_b, "color", "#295D8A")
    setTrackStyleParam(sgRNA_JQ1_FCPF_high_b, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_FCPF_high_b, "ylabgp", list("cex"=0.6, fontface="bold"))
    
    sgRNA_JQ1_high <- importScore(bwfiles["5a"],  format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_high, "color", "#9400D3")
    setTrackStyleParam(sgRNA_JQ1_high, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_high, "ylabgp", list("cex"=0.6, fontface="bold"))

    sgRNA_JQ1_high_b <- importScore(bwfiles["5b"],  format="BigWig", ranges=gr)
    setTrackStyleParam(sgRNA_JQ1_high_b, "color", "#9400D3")
    setTrackStyleParam(sgRNA_JQ1_high_b, "height", .1998)
    setTrackStyleParam(sgRNA_JQ1_high_b, "ylabgp", list("cex"=0.6, fontface="bold"))


    tl <- trackList(
        gTrack,
                    #sgRNA_JQ1_high_b,
                    sgRNA_JQ1_high,
                    #sgRNA_JQ1_FCPF_high_b,
                    sgRNA_JQ1_FCPF_high,
                    #sgRNA_JQ1_FCPF_low_b,
                    sgRNA_JQ1_FCPF_low,
                    #JQ1_FCPF_low_b,
                    JQ1_FCPF_low,
                    NTb#,
                   # NT
    )
    
    names(tl) <- c(
        symbol,
                    #"sgRNA-JQ1 high_b",
                    "sgRNA-JQ1 high",
                    #"sgRNA-JQ1_FCPF high_b",
                    "sgRNA-JQ1_FCPF high",
                    #"sgRNA-JQ1_FCPF low_b",
                    "sgRNA-JQ1_FCPF low",
                    #"JQ1-FCPF low_b",
                    "JQ1-FCPF low",
                    "NT "#b",
                   # "NT"
    )

    # get max score
    max <- 0
    for(i in 1:length(tl)){
        #print(tl[[i]]$dat$score)
        new <- max(tl[[i]]$dat$score)
        if (new > max){
            max <- new
        }
        if (new > max_overall){
            max_overall <- new
        }
        
    }


    # fixed y-axis
    for(j in 2:length(tl)){
        setTrackStyleParam(tl[[j]], "ylim", c(0, max))
    }
    # flip y-axis position
    # for(i in 1:length(tl)){
    #     setTrackYaxisParam(tl[[i]], "main", FALSE)
    # #         setTrackYaxisParam(NT, "label", T)
    # # setTrackYaxisParam(NT, "at", c(0.2, 0.8))
    # # setTrackYaxisParam(NT, "gp", list("cex"=0.6))

    # }

    # svg(paste0("output/tracks/", symbol ,".svg"), width= 3000, height = 2200)
    # viewTracks(tl, gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)
    # dev.off()
    png(paste0("output/tracks/", symbol ,"_special.png"), width=3000, height = 2200, res=300)
    viewTracks(tl, gr=gr, viewerStyle=viewerStyle) #, smooth=T)
    dev.off()

    svg(paste0("output/tracks/", symbol ,"_special.svg"), width=10, height = 7.3)
    viewTracks(tl, gr=gr, viewerStyle=viewerStyle) #, smooth=T)
    dev.off()
   

}

