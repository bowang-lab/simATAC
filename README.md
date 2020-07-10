# simATAC

simATAC is an R package for the simple simulation of single-cell ATAC sequencing 
(scATAC-seq) data. simATAC simulates feature matrices from real scATAC-seq bin by cell 
matrix by estimating the modelling parameters of three main metrics: library size, 
bin read count mean, and non-zero cells proportion. simATAC deploys the estimated 
parameters to simualte final synthetic counts. 

simATAC stores simulation parameters in a simATACCount class. It is built on top of 
[`scater`][scater] and stores simulated values in [`SingleCellExperiment`][SCE] objects. 

The [`SingleCellExperiment`][SCE] object generated by simATAC can be converted to 
binary count matrix, peak by cell matrix, and any feature matrix version by giving 
desired genomic regions in BED file.

## Installation.


## Getting started

We provided simATAC tutorial with examples and explanations about its functions and how to use them. This is a detailed document that introduces the main features of simATAC.



[scater]: https://github.com/davismcc/scater
[SCE]: https://github.com/drisso/SingleCellExperiment
[contrib]: https://github.com/Bioconductor/Contributions/issues/209
[bioc]: https://bioconductor.org/packages/devel/bioc/html/splatter.html
[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
[paper]: http://dx.doi.org/10.1186/s13059-017-1305-0
