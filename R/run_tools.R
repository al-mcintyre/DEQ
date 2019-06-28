run.tools <- function(results,peak.counts,meta,tool,input.bams,ip.bams,treated.input.bams,treated.ip.bams){
  tools <- strsplit(tool,'')[[1]]
  full.results <- TRUE
  if ('d' %in% tools){
    results <- cbind(results,run.deseq2(peak.counts,meta))
  }
  if ('e' %in% tools){
    results <- cbind(results,run.edger(peak.counts,meta))
  }
  if ('q' %in% tools){
    results <- cbind(results,run.qnb(peak.counts,meta))
  }
  if ('m' %in% tools){
    results <- cbind(results,run.metdiff(peak.counts,meta))
  }
  return(results)
}

run.deseq2 <- function(cnts,meta){
  inf.dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnts,colData = meta,design = ~Condition+IP+Condition:IP)
  inf.dds.LRT <- DESeq2::DESeq(inf.dds,betaPrior=FALSE, test="LRT",
                       full=~Condition+IP+Condition:IP,reduced=~Condition+IP)    
  inf.dds.res <- DESeq2::results(inf.dds.LRT)
  results <- as.data.frame(cbind(inf.dds.res$pvalue,inf.dds.res$padj))
  colnames(results) <- c("deseq2.p","deseq2.padj")
  return(results)
}

run.edger <- function(cnts,meta){
  #add count filter?
  er.design <- model.matrix(~meta$Condition+meta$IP+meta$Condition*meta$IP)
  er.dgelist <- edgeR::DGEList(counts=cnts,group=meta$Condition) 
  er.dgelist <- edgeR::estimateDisp(er.dgelist, design=er.design)
  er.fit <- edgeR::glmFit(er.dgelist, er.design)
  er.lrt <- edgeR::glmLRT(er.fit, coef=4)
  #hist(er.lrt$table$PValue) er.lrt$table$logFC,
  results <- as.data.frame(cbind(er.lrt$table$PValue,p.adjust(er.lrt$table$PValue)))
  colnames(results) <- c("edger.p","edger.padj")
  return(results)
}

run.qnb <- function(cnts,meta){
  meth1 <- cnts[,which(meta$Condition == 'treatment' & meta$IP == "IP")]
  meth2 <- cnts[,which(meta$Condition != 'treatment' & meta$IP == "IP")]
  unmeth1 <- cnts[,which(meta$Condition == 'treatment' & meta$IP == "input")]
  unmeth2 <- cnts[,which(meta$Condition != 'treatment' & meta$IP == "input")]
  qnb.result.sim = QNB::qnbtest(meth1, meth2, unmeth1, unmeth2, mode="per-condition")
  #qnb.result.sim$log2.RR,qnb.result.sim$log2.OR,
  results <- as.data.frame(cbind(qnb.result.sim$pvalue,qnb.result.sim$padj))
  colnames(results) <- c("qnb.p","qnb.padj")
  return(results)
}

run.metdiff <- function(cnts,meta){
  meth1 <- cnts[,which(meta$Condition == 'treatment' & meta$IP == "IP")]
  meth2 <- cnts[,which(meta$Condition != 'treatment' & meta$IP == "IP")]
  unmeth1 <- cnts[,which(meta$Condition == 'treatment' & meta$IP == "input")]
  unmeth2 <- cnts[,which(meta$Condition != 'treatment' & meta$IP == "input")]
  metdiff.result <- diff.call.module(meth1,unmeth1,meth2,unmeth2)
  results <- as.data.frame(cbind(metdiff.result$DIFF$pvalues,metdiff.result$DIFF$fdr))
  colnames(results) <- c("metdiff.p","metdiff.padj")
  return(results)
}

#for genes l2FC and peak IP l2FC
run.deseq2.4l2fc <- function(cnts,meta,label){
  dds <- DESeq2::DESeqDataSetFromMatrix(cnts,meta,formula(~Condition))
  dds$Condition <- factor(dds$Condition, levels=c('control','treatment'))
  gene.col2check <- meta$Condition
  dds$Condition <- droplevels(dds$Condition)
  gene.deseq <- DESeq2::DESeq(dds)
  gene.deseq <- DESeq2::results(gene.deseq)
  gene.results <- gene.deseq[,c("log2FoldChange","pvalue","padj")]
  colnames(gene.results) <- paste0(label,c(".l2fc",".p",".padj"))
  return(gene.results)
}
