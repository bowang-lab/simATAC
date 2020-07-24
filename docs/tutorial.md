Welcome to simATAC simulator!
simATAC is an R package developed for simualting single-cell ATAC sequencing (scATAC-seq) count matrices. simATAC is an easy to use simulation framework that given a group of cells having similar biological characteristics in the format of bin by cell matrix as an input, it generates a synthetic bin by cell matrix resembling the real samples. simATAC mainly performs two estimation and simualtion steps in order to generate final synthetic counts. There are separate functions for each step, and this tutorial gives an overview and introduction to simATAC’s functionality.

Assuming there is scATAC-seq dataset with a sepcific cell type, you can convert raw BAM files into a bin by cell (5kbp-window) matrix with any customized pipeline. We used Snaptools to generate .snap files which contains bin by cell array in it. In order to run examples, we use the snap file of GSE99172 real scATAC-seq sample. We will skip the snap generation (See [here](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap) for how to generate a snap file). Instead, we will download the snap file, which is provided in the example folder. 


## Table of Contents

- [Loading simATAC](#load)
- [simATACCount class](#simATACCount)
- [Estimation function](#estimation)
- [Simulation function](#simulation)
- [Set parameters](#set)
- [Get parameters](#get)


<a name="simATACCount"></a>**Loading simATAC**        
We need to load simATAC R package to be able to use the provided functions.

```bash
$ library(simATAC)
```

<a name="simATACCount"></a>**simATACCount class**        
We defined a simATACCount class, which is an class specifically desinged for storing simATAC scATAC-seq simulation parameters. You can create a new simATACCount object by:

```bash
$ object <- newsimATACCount()
```
and print its parameters:


