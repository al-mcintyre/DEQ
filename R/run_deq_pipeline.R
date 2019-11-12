#' DEQ 
#' 
#' A package to conveniently run DESeq2, edgeR, and QNB for the detection of differential methylation in MeRIP/m6A-seq data.
#'
#' @param input.bams vector of control input bam files
#' @param ip.bams vector of control IP bam files
#' @param input.bams vector of treatment input bam files
#' @param ip.bams vector of treatment IP bam files
#' @param peak.files vector of peak file names
#' @param gtf location of a gtf file for annotations
#' @param paired.end boolean indicating whether reads were paired end (defaults
#'   to FALSE)
#' @param outfi output file name (defaults to deq.results.txt)
#' @param tool which tools to run (options: any combination of d, e, and q;
#'   defaults to all: 'deq')
#' @param compare.gene whether to calculate changes in gene expression for
#'   comparison (highly recommended, defaults to TRUE)
#' @param readlen average read length for the data (defaults to 100)
#' @param fraglen average fragment length for the data (defaults to 100) --
#'   (fraglen-readlen) is used to calculate the extension of reads for counting
#' @param nthreads number of threads to run on (defaults to 1)
#' @return GRanges object with peaks and estimated changes
#' @import dplyr 
#' @import GenomeInfoDb
#' @import S4Vectors
#' @export
#' @include import_peaks.R count_reads.R run_tools.R metdiff_function_copy.R
#' @examples 

deq <- function(input.bams,ip.bams,treated.input.bams,treated.ip.bams,
                peak.files,gtf,paired.end=FALSE,outfi='deq_results.txt',
                tool='deq',compare.gene=TRUE,readlen=100,fraglen=100,nthreads=1){
  
  if (length(input.bams) != length(ip.bams) | length(treated.input.bams) != length(treated.ip.bams)){
    stop('number of IP bam files must equal number of input bam files')
  }
  n.c <- length(input.bams)
  n.t <- length(treated.input.bams)
  extension <- fraglen-readlen
  meta.data <- data.frame(Condition=c(rep('control',n.c*2),rep('treatment',n.t*2)),
                     IP=c(rep('input',n.c),rep('IP',n.c),rep('input',n.t),rep('IP',n.t)))
  
  #load gtf annotations
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf,format='gtf')
  gtf.in <- rtracklayer::import(gtf)
  genenames <- unique(as.data.frame(gtf.in)[,c("gene_id","gene_name","strand")])
  anno <- list(txdb=txdb,genenames=genenames)
  annotation.order <- c("utr3","utr5","exon","intron","utr3*","utr5*")
  
  #load peaks and annotate
  peaks <- import.peaks(peak.files,anno,annotation.order) 
  
  #count reads
  all.bams <- c(input.bams,ip.bams,treated.input.bams,treated.ip.bams)
  peaks <- count.reads(peaks,all.bams,paired.end,extension)
  peak.counts <- DESeq2::counts(peaks$peak.counts)
  
  #run DESeq2, edgeR, and QNB to predict changes in m6A methylation
  results <- peaks$peaks[,c("annot","main.gene")]
  results <- run.tools(results,peak.counts,meta.data,tool,input.bams,ip.bams,treated.input.bams,treated.ip.bams) 
  peaks$peak.de <- run.deseq2.4l2fc(peak.counts[,which(meta.data$IP == "IP")],
                                  meta.data[which(meta.data$IP == "IP"),],'peak')
  results$peak.l2fc <- peaks$peak.de$peak.l2fc
  
  #calculate gene log2 fold change  
  if (compare.gene){
    peaks$gene.counts <- get.gene.counts(c(input.bams,treated.input.bams),
                                         gtf,paired.end,extension,genenames)
    peaks$gene.de <- run.deseq2.4l2fc(peaks$gene.counts,meta.data[which(meta.data$IP == "input"),],'gene')
    results$gene.l2fc <- peaks$gene.de[results$main.gene,]$gene.l2fc
    results$gene.p <- peaks$gene.de[results$main.gene,]$gene.p
    results$gene.padj <- peaks$gene.de[results$main.gene,]$gene.padj
    results$diff.l2fc <- results$peak.l2fc - results$gene.l2fc 
  }
  
  results$start <- results$start-1
  colnames(peak.counts) <- make.names(paste0(meta.data$Condition,"_",meta.data$IP),unique=TRUE)
  write.table(peak.counts,gsub('.txt','.counts.txt',outfi),quote = FALSE,sep='\t',row.names = TRUE,col.names=TRUE)
  write.table(results,outfi,quote = FALSE,sep='\t',row.names = FALSE,col.names=TRUE)
  return(results)
}
