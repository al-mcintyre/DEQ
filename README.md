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

This package has been tested on MacOS (v10.13.6) and Linux (Red Hat Enterprise Linux 6.3), using R v3.6.0 and v3.6.1. Typical install time < 1 min

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

Demo data is provided. After downloading, run
``` 
input.bams <- c('untreated_total_input_1.chr21.star.sorted.bam','untreated_total_input_2.chr21.star.sorted.bam')
ip.bams <- c('untreated_total_IP_1.chr21.star.sorted.bam','untreated_total_IP_2.chr21.star.sorted.bam')
treated.input.bams <- c('heatshock_total_input_1.chr21.star.sorted.bam','heatshock_total_input_2.chr21.star.sorted.bam')
treated.ip.bams <- c('heatshock_total_IP_1.chr21.star.sorted.bam','heatshock_total_IP_2.chr21.star.sorted.bam')
peak.files <- 'peaks.bed'
#gtf <- location of hg38 gtf file
deq(input.bams, ip.bams, treated.input.bams, treated.ip.bams, peak.files,
  gtf, paired.end = FALSE, outfi = "deq_results.txt", tool = "deq",
  compare.gene = TRUE, readlen = 50, fraglen = 100, nthreads = 1)
  ```
from the appropriate folder (or adjust file paths). Must provide a gtf file for hg38. Run-time should be ~5 minutes for this demo. 

Further scripts to generate figures for our [bioRxiv paper](https://www.biorxiv.org/content/10.1101/657130v2) are available [here](https://github.com/al-mcintyre/merip_reanalysis_scripts).
