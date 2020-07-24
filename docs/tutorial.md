Welcome to simATAC simulator!
simATAC is an R package developed for simualting single-cell ATAC sequencing (scATAC-seq) count matrices. simATAC is an easy to use simulation framework that given a group of cells having similar biological characteristics in the format of bin by cell matrix as an input, it generates a synthetic bin by cell matrix. simATAC mainly performs two estimation and simualtion steps in order to generate final synthetic counts. There are separate functions for each step, and this tutorial gives an overview and introduction to simATAC’s functionality.

Assuming there is scATAC-seq dataset with a sepcific cell type, you can convert raw BAM files into a bin by cell matrix with any pipeline. We used Snaptools to generate .snap files which contains bin by cell array. In order to run examples, we use the snap file of GSE99172 real scATAC-seq sample, which is provided in the example folder. 


## Table of Contents

- [simATACCount class](#simATACCount)
- [Estimation step](#estimation)
- [Simulation step](#simulation)
- [Set parameters](#set)
- [Get parameters](#get)


<a name="simATACCount"></a>**simATACCount class**        
In this example, we will skip the snap generation (See [here](https://github.com/r3fang/SnapATAC/blob/master/examples/10X_PBMC_15K/README.md) for how to generate a snap file). Instead, we will download the snap file. The downloaded snap file already contains the cell-by-bin/cell-by-peak matrix.

```bash
$ wget http://renlab.sdsc.edu/r3fang/share/github/Mouse_Brain_10X/atac_v1_adult_brain_fresh_5k.snap
$ http://renlab.sdsc.edu/r3fang/share/github/Mouse_Brain_10X/atac_v1_adult_brain_fresh_5k_singlecell.csv
```
