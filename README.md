# DEQ
DEteQtion of changes in m6a. An R package to conveniently run DESeq2, edgeR, and QNB for the detection of differential methylation in MeRIP/m6A-seq data.

### Installation

Install in R using `devtools::install_github("al-mcintyre/DEQ")`

May need to import S4Vectors separately using `library(S4Vectors)`

Need to install dependencies separately for now. 
```
BiocManager::install("DESeq2")
BiocManager::install("edgeR")
BiocManager::install("GenomicFeatures")
BiocManager::install("Rsamtools")
```
QNB is no longer supported on CRAN, and requires exomePeak. It can be installed from source. 
```
BiocManager::install("exomePeak")
install.packages("https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz", repos = NULL, type="source")
```

This package has been tested on MacOS (v10.13.6) and Linux (Red Hat Enterprise Linux 6.3), using R v3.6.0 and v3.6.1

### Running DEQ

Run using 
``` 
deq(input.bams, ip.bams, treated.input.bams, treated.ip.bams, peak.files,
  gtf, paired.end = FALSE, outfi = "deq_results.txt", tool = "deq",
  compare.gene = TRUE, readlen = 100, fraglen = 100, nthreads = 1)
  ```
(see ?deq for details)

Output:
1. Counts file: number of reads per peak for each input bam file (labelled as control or treatment replicates)
2. Results file: chromosomal locations and annotations for peaks, along with significance predictions for any tools run (with and without p value adjustment), and log2 fold changes for IP reads within peaks (peak.l2fc) and input reads for associated genes (gene.l2fc) calculated using DESeq2, significance for gene expression changes between conditions from DESeq2, and the difference between peak.l2fc and gene.l2fc (diff.l2fc)
