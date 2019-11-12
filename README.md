# DEQ
DEteQtion of changes in m6a. An R package to conveniently run DESeq2, edgeR, and QNB for the detection of differential methylation in MeRIP/m6A-seq data.

Install in R using `devtools::install_github("al-mcintyre/DEQ")`

Run using 
``` 
deq(input.bams, ip.bams, treated.input.bams, treated.ip.bams, peak.files,
  gtf, paired.end = FALSE, outfi = "deq_results.txt", tool = "deq",
  compare.gene = TRUE, readlen = 100, fraglen = 100, nthreads = 1)
  ```
(see ?deq for details)

May need to import S4Vectors separately - I've had some issues `library(S4Vectors)`

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
