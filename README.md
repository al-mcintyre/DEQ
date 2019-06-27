# DEQ
An R package to conveniently run DESeq2, edgeR, and QNB for the detection of differential methylation in MeRIP/m6A-seq data.

Install in R using `devtools::install_github("al-mcintyre/DEQ")`

Run using 
``` 
deq(input.bams, ip.bams, treated.input.bams, treated.ip.bams, peak.files,
  gtf, paired.end = FALSE, outfi = "deq_results.txt", tool = "deq",
  compare.gene = TRUE, readlen = 100, fraglen = 100, nthreads = 1)
  ```
(see ?deq for details)
