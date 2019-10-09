get.gene.counts <- function(bamlist,gtffi,pe,extension,genenames){
  gene.data <- Rsubread::featureCounts(files = bamlist, GTF.featureType="exon",GTF.attrType="gene_id",
                                         annot.ext=gtffi, isGTFAnnotationFile=TRUE,isPairedEnd = pe,readExtension3 = extension)
  gene.counts <- gene.data$counts
  conversion.table <- as.data.frame(genenames)
  row.names(conversion.table) <- make.names(genenames$gene_id,unique=TRUE)
  row.names(gene.counts) <- conversion.table[row.names(gene.counts),"gene_name"]
  return(gene.counts)
}

#following two functions also adapted from CLIPanalyze
checkBamFileList <- function (bamfiles, clean.names = TRUE) {
  if (is.character(bamfiles)) {
    bamfiles <- Rsamtools::BamFileList(bamfiles)
  }
  if (!is(bamfiles, "BamFileList")) {
    stop("BamFileList required")
  }
  if (!all(file.exists(BiocGenerics::path(bamfiles)))) {
    lost <- !file.exists(BiocGenerics::path(bamfiles))
    stop(paste(BiocGenerics::path(bamfiles[lost]), collapse = ","), " BAMs not found")
  }
  if (clean.names) {
    nms <- names(bamfiles)
    if (is.null(nms)) {
      nms <- BiocGenerics::path(bamfiles)
    }
    nms <- gsub("\\.bam$", "", basename(nms))
    names(bamfiles) <- make.unique(nms)
  }
  bamfiles
}

count.reads <- function(peak.data,bamfiles,paired.end,extension,sample.names = NULL,
                        stranded = 0,nthreads=1){
    message("count reads in peaks...")
    bamfiles.filenames <- bamfiles
    bamfiles <- checkBamFileList(bamfiles)
    if (is.list(peak.data)) {
      peaks <- peak.data$peaks
      if (!("peak.counts" %in% names(peak.data))) {
        peak.data$peak.counts <- NULL
      }
      if (!("gene.counts.nopeaks" %in% names(peak.data))) {
        peak.data$gene.counts.nopeaks <- NULL
      }
    }
    else {
      peaks <- peak.data
      peak.data <- list(peaks = peaks, peak.counts = NULL, 
                        gene.counts.nopeaks = NULL)
    }
    if (nthreads > 1) {
      message(sprintf("use %s threads in parallel", nthreads))
      BiocParallel::register(BiocParallel::MulticoreParam(workers = nthreads))
    }
    peaks.for.counts <- peaks
    peak.df <- data.frame(GeneID = names(peaks.for.counts), Chr = as.character(GenomicRanges::seqnames(peaks.for.counts)), 
                          Start = GenomicRanges::start(peaks.for.counts), End = GenomicRanges::end(peaks.for.counts), 
                          Strand = GenomicRanges::strand(peaks.for.counts))
    peak.counts <- Rsubread::featureCounts(bamfiles.filenames, 
                                           annot.ext = peak.df, allowMultiOverlap = TRUE, 
                                           nthreads = nthreads, minFragLength = 20, isPairedEnd = paired.end, readExtension3 = extension)
    #strandSpecific = stranded, 
    peak.counts <- peak.counts$counts
    if (!is.null(sample.names)) {
      if (length(colnames(peak.counts)) == length(sample.names)) {
        colnames(peak.counts) <- sample.names
      }
      else {
        warning("length of sample.names does not match colnames(peak.counts)")
      }
    }
    peak.counts <- DESeq2::DESeqDataSetFromMatrix(peak.counts, 
                                                  design = as.formula("~1"), colData = data.frame(row.names = colnames(peak.counts)))
    peak.data$peak.counts <- peak.counts

    return(peak.data)
}
