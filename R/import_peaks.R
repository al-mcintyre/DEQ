import.peaks <- function(peak.files,anno,annot.order=annotation.order){
  if (exists("loaded.peaks")){ remove(loaded.peaks)}
  for (file in peak.files){
    if (exists("loaded.peaks")){
      loaded.peaks <- GenomicRanges::union(loaded.peaks,rtracklayer::import(file,format="BED"))
    }else{
      loaded.peaks <- rtracklayer::import(file,format="BED")
    }
  }
  loaded.peaks <- GenomicRanges::reduce(loaded.peaks) #merge overlapping peaks
  loaded.peaks <- GenomeInfoDb::keepStandardChromosomes(loaded.peaks,pruning.mode = "coarse")
  #annotate macs2 peaks
  peaks.anno <- annotate.peaks(loaded.peaks,anno,annot.order)
  remove(loaded.peaks)
  return(peaks.anno)
}

annotate.peaks <- function(loaded.peaks,anno,annot.order){
  names(loaded.peaks) <- c(paste0("peak",1:length(loaded.peaks)))
  peaks.anno <- annotateFeatures(loaded.peaks, txdb = anno$txdb,
                                   genenames = anno$genenames,
                                   utr5.extend = 2000,
                                   utr3.extend = 2000)
  names(peaks.anno) <- c("peaks", "overlaps")
  peaks.anno$peaks <- prioritizeAnnotations(peaks.anno$peaks,feature.order = annot.order, genenames = anno$genenames) 
  return(peaks.anno)
}

#the following functions were adapted from CLIPanalyze
findOverlapsFeature <- function (gr, features, gene.ids, feature.name){
  hits <- GenomicRanges::findOverlaps(gr, features)
  hits <- data.table::as.data.table(hits)
  hits$tx_name <- names(features)[hits$subjectHits]
  hits <- merge(hits, gene.ids, by = "tx_name")
  hits <- hits[, c("queryHits", "gene_name")] %>% unique(by = NULL)
  hits.arbitrary <- hits %>% unique(by = "queryHits")
  hits$name <- gr[hits$queryHits]$name
  hits <- hits[, c("name", "gene_name")] %>% unique(by = NULL)
  mcols(gr)[, feature.name] <- as.character(NA)
  mcols(gr)[hits.arbitrary$queryHits, feature.name] <- hits.arbitrary$gene_name
  return(list(gr = gr, hits = hits))
}

annotateFeatures <- function (gr, txdb, genenames, utr5.extend = 2000, utr3.extend = 2000){
  gr$name <- names(gr)
  message("prepare transcript and gene names... ")
  tx.names <- AnnotationDbi::keys(txdb, "TXNAME")
  gene.ids <- AnnotationDbi::select(txdb, keys = tx.names, 
                                    columns = "GENEID", keytype = "TXNAME")
  colnames(gene.ids) <- c("tx_name", "gene_id")
  gene.ids <- data.table::as.data.table(gene.ids)
  genenames <- data.table::as.data.table(genenames)
  gene.ids <- merge(gene.ids, genenames, by = "gene_id")
  all.overlaps <- list()
  message("annotate with exons... ")
  exons <- GenomicFeatures::exonsBy(txdb, use.names = TRUE)
  gr.hits <- findOverlapsFeature(gr, exons, gene.ids, "exon")
  gr <- gr.hits[["gr"]]
  all.overlaps <- c(all.overlaps, list(exon = gr.hits[["hits"]]))
  message("annotate with introns... ")
  introns <- GenomicFeatures::intronsByTranscript(txdb, use.names = TRUE)
  gr.hits <- findOverlapsFeature(gr, introns, gene.ids, "intron")
  gr <- gr.hits[["gr"]]
  all.overlaps <- c(all.overlaps, list(intron = gr.hits[["hits"]]))
  message("annotate with 5'UTRs... ")
  five.utr <- GenomicFeatures::fiveUTRsByTranscript(txdb, use.names = TRUE)
  gr.hits <- findOverlapsFeature(gr, five.utr, gene.ids, "utr5")
  gr <- gr.hits[["gr"]]
  all.overlaps <- c(all.overlaps, list(utr5 = gr.hits[["hits"]]))
  if (!is.null(utr5.extend) && !is.na(utr5.extend) && utr5.extend > 
      0) {
    message(sprintf("annotate as utr5* if %s upstream of 5'UTRs... ", 
                    utr5.extend))
    five.utr.upstream <- GenomicRanges::flank(five.utr, width = utr5.extend, 
                                              start = TRUE)
    gr.hits <- findOverlapsFeature(gr, five.utr.upstream, 
                                   gene.ids, "utr5*")
    gr <- gr.hits[["gr"]]
    all.overlaps <- c(all.overlaps, list(`utr5*` = gr.hits[["hits"]]))
  }
  message("annotate with 3'UTRs... ")
  three.utr <- GenomicFeatures::threeUTRsByTranscript(txdb, 
                                                      use.names = TRUE)
  gr.hits <- findOverlapsFeature(gr, three.utr, gene.ids, "utr3")
  gr <- gr.hits[["gr"]]
  all.overlaps <- c(all.overlaps, list(utr3 = gr.hits[["hits"]]))
  if (!is.null(utr3.extend) && !is.na(utr3.extend) && utr3.extend > 
      0) {
    message(sprintf("annotate as utr3* if %s downstream of 3'UTRs... ", 
                    utr3.extend))
    three.utr.upstream <- GenomicRanges::flank(three.utr, 
                                               width = utr3.extend, start = FALSE)
    gr.hits <- findOverlapsFeature(gr, three.utr.upstream, 
                                   gene.ids, "utr3*")
    gr <- gr.hits[["gr"]]
    all.overlaps <- c(all.overlaps, list(`utr3*` = gr.hits[["hits"]]))
  }
  message("annotate with intergenic... ")
  if (length(gr) > 0) {
    gr$intergenic <- TRUE
    which.intergenic <- !is.na(gr$exon) | !is.na(gr$intron)
    if (any(which.intergenic)) {
      gr[which.intergenic]$intergenic <- FALSE
    }
    if (("utr5*" %in% colnames(mcols(gr))) && any(!is.na(gr$"utr5*"))) {
      gr[!is.na(gr$"utr5*")]$intergenic <- FALSE
    }
    if (("utr3*" %in% colnames(mcols(gr))) && any(!is.na(gr$"utr3*"))) {
      gr[!is.na(gr$"utr3*")]$intergenic <- FALSE
    }
  }
  all.overlaps <- plyr::ldply(all.overlaps, .id = "feature")
  all.overlaps <- data.table::as.data.table(all.overlaps)
  return(list(gr = gr, overlaps = all.overlaps))
}

prioritizeAnnotations <- function (gr, feature.order, genenames){
  message(sprintf("prioritize annotations in this order: %s ...", 
                  paste0(feature.order, collapse = ", ")))
  feature.diff <- setdiff(feature.order, colnames(mcols(gr)))
  if (length(feature.diff) > 0) {
    message(sprintf("warning: these annotations are not avaiable: %s", 
                    paste0(feature.diff, collapse = ", ")))
  }
  feature.order <- feature.order[feature.order %in% colnames(mcols(gr))]
  mcols(gr)[, "annot"] <- as.character(NA)
  for (feature in feature.order) {
    assign.feature <- is.na(mcols(gr)[, "annot"]) & !is.na(mcols(gr)[, 
                                                                     feature])
    mcols(gr)[assign.feature, "annot"] <- feature
    mcols(gr)[assign.feature, "main.gene"] <- mcols(gr)[assign.feature,feature]
  }
  # get strands based on primary gene annotations
  strand.table <- as.data.frame(genenames)
  row.names(strand.table) <- make.names(genenames$gene_name,unique=TRUE)
  assign.feature <- !is.na(mcols(gr)[, "annot"]) & mcols(gr)[, "annot"] != "intergenic"
  strand(gr[assign.feature,]) <- strand.table[gsub('-','.',mcols(gr)[assign.feature,"main.gene"]),"strand"]
  mcols(gr)[assign.feature, "main.gene.id"] <- strand.table[gsub('-','.',mcols(gr)[assign.feature,"main.gene"]),"gene_id"]
  
  assign.feature <- is.na(mcols(gr)[, "annot"])
  mcols(gr)[assign.feature, "annot"] <- "intergenic"
  strand(gr[assign.feature,]) <- '+' #necessary approximation for downstream functions
  return(gr)
}
