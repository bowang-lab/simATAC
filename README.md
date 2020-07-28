# simATAC

simATAC is a framework provided as an R package that generates a single-cell Assay for Transposase-Accessible Chromatin sequencing (scATAC-seq) count matrix, highly resembling real scATAC-seq datasets in library size, sparsity, and averaged chromatin accessibility signals. simATAC deploys statistical functions derived from analyzing 90 real scATAC-seq cell groups to model read distributions. simATAC provides a robust and systematic approach to generate \textit{in silico} scATAC-seq samples with cell labels for a comprehensive tool assessment.


## Installation.

Install simATAC R package (R >= 3.5.0): 

```bash
$ R
> library(devtools)
> install_github("bowang-lab/simATAC")
```

## Getting started


simATAC tutorial is provided with examples and explanations about its functions and how to use them. This is detailed documentation that introduces the main features of simATAC.
* [simATAC tutorial](https://github.com/bowang-lab/simATAC/blob/master/docs/tutorial.md)




[scater]: https://github.com/davismcc/scater
[SCE]: https://github.com/drisso/SingleCellExperiment
[contrib]: https://github.com/Bioconductor/Contributions/issues/209
[bioc]: https://bioconductor.org/packages/devel/bioc/html/splatter.html
[vignette]: https://bioconductor.org/packages/devel/bioc/vignettes/splatter/inst/doc/splatter.html
[paper]: http://dx.doi.org/10.1186/s13059-017-1305-0
