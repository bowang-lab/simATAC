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


<a name="load"></a>**Loading simATAC**        
We need to load simATAC R package to be able to use the provided functions.

```bash
> library(simATAC)
```

<a name="simATACCount"></a>**simATACCount class**        
We defined a simATACCount class, which is an class specifically desinged for storing simATAC scATAC-seq simulation parameters. You can create a new simATACCount object by:

```bash
> object <- newsimATACCount()
> object
```

which will print the default values for its parameters:

```bash
An object of class "simATACCount"
Slot "nBins":
[1] 642098

Slot "nCells":
[1] 500

Slot "seed":
[1] 425649

Slot "default":
[1] TRUE

Slot "species":
[1] "human"

Slot "lib.mean1":
[1] 13.60503

Slot "lib.mean2":
[1] 14.93826

Slot "lib.sd1":
[1] 1.745264

Slot "lib.sd2":
[1] 1.009923

Slot "lib.prob":
[1] 0.5257138

Slot "non.zero.pro":
[1] 1

Slot "mean.coef0":
[1] 0.002822035

Slot "mean.coef1":
[1] 0.6218985

Slot "mean.coef2":
[1] 1.976122

Slot "noise.mean":
[1] 0

Slot "noise.sd":
[1] 0
and print the description for its parameters:

```

You can print the description for each parameter documented in the R package by running

```bash
> ?simATACCount
```

```bash
simATACCount              package:simATAC              R Documentation

The simATACCount class

Description:

     S4 class that holds parameters for the count matrix of simATAC
     simulation.

Parameters:

     simATAC simulation parameters:

     ‘nBins’ The bin number to simulate.

     ‘nCells’ The cell number to simulate.

     ‘[seed]’ Seed to use for generating random numbers.

     ‘[default]’ The logical variable whether to use default parameters
          (TRUE) or learn from data (FALSE)

     ‘[species]’ The string indicating the species of the input cells

     _Library size parameters_

          ‘lib.mean1’ Mean parameter for the first component of library
              size bimodal Gaussian distribution.

          ‘lib.mean2’ Mean parameter for the second component of
              library size bimodal Gaussian distribution

          ‘lib.sd1’ Standard deviation parameter for the first
              component of library size bimodal Gaussian distribution.

          ‘lib.sd2’ Standard deviation parameter for the second
              component of library size bimodal Gaussian distribution.

          ‘lib.prob’ Probability parameter for the first component in
              bimodal Gaussian distribution. The probability for the
              second component is 1-lib.prob.

     _Zero entry parameters_

          ‘non.zero.pro’ The proportion of non-zero cells per bin in
              the original count matrix

          ‘mean.coef0’ Estimated coefficient of power zero variable in
              the polynomial function

          ‘mean.coef1’ Estimated coefficient of power one variable in
              the polynomial function

          ‘mean.coef2’ Estimated coefficient of power two variable in
              the polynomial function

     _[noise]_

          ‘[noise.mean]’ Gaussian mean to be added as noise to the
              final simulated counts

          ‘[noise.sd]’ Gaussian standard deviation to be added as noise
              to the final simulated counts

     The parameters not shown in brackets can be estimated from real
     data using ‘simATACEstimate’. For details of the simATAC
     simulation see ‘simATACSimulate’.
```


<a name="estimation"></a>**Estimation function**


simATAC generates a synthetic scATAC-seq count matrix by first fitting statistical models to the three main parameters, including library size (read coverage of cells), bin mean (the average of counts per bin), and bin non-zero cell proportion (non-zero cell proportion in each bin). simATAC builds upon Gaussian mixture distribution to model cell library sizes, and polynomial regression model to represent the relationship between the bin means and the non-zero cell proportions of bins. 
simATAC allows you to estimate the parameters of the parameters' models by simATACEstimate() function:

```bash
## return the cell by bin matrix from the snap file
> count <- getCountFromh5("GSE99172.snap")
## print dimentionality of the count matrix
> dim(count)
[1]    288 642098
```
count object includes 288 cells, and 642098 bins with 5000 base pair length. It is a sparse matrix abnd you can print the type of it by running
```bash
> typeof(count)
[1] "S4"
> class(count)
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"
```

Because the ocunt object is cell by bin (rows are cells and columns are bins), you need to convert it to the bin by cell to be able to feed it into simATAC package.

```bash
library(Matrix)

> object <- simATACEstimate(t(count))
simATAC is:
...Estimating library size...
...Estimating non-zero cell proportion...
...Estimating bin mean...
>
> object
An object of class "simATACCount"
Slot "nBins":
[1] 642098

Slot "nCells":
[1] 500

Slot "seed":
[1] 425649

Slot "default":
[1] FALSE

Slot "species":
[1] "human"

Slot "lib.mean1":
[1] 13.6052

Slot "lib.mean2":
[1] 14.93831

Slot "lib.sd1":
[1] 1.745244

Slot "lib.sd2":
[1] 1.009877

Slot "lib.prob":
[1] 0.5257947

Slot "non.zero.pro":
    [1] 0.000000000 0.000000000 0.003472222 0.006944444 0.000000000 0.003472222
    [7] 0.003472222 0.000000000 0.000000000 0.000000000 0.006944444 0.003472222
   [13] 0.000000000 0.000000000 0.003472222 0.006944444 0.000000000 0.000000000
    [ reached getOption("max.print") -- omitted 642080 entries ]
    
Slot "mean.coef0":
[1] 0.002822035

Slot "mean.coef1":
[1] 0.6218985

Slot "mean.coef2":
[1] 1.976122

Slot "noise.mean":
[1] 0

Slot "noise.sd":
[1] 0
```

simATACEstimate function estimates the paramters of fitted models and if the verbose vairable is set to TRUE (which is by default), it prints the progress of estimation process. 
1. simATAC first fits a Guassian mixture model to the log-transformed of library size.
1. It then calculates the proportion of cell having non-zero count whitin each bin (non-zero cell proportion).
1. Finally fits a polynomial regression model to the relation between bin means and bin non zero cell proportions and estiamtes the parameters of the polynomial function.

All estimated model parameters are stored in the simATACCount object. 

The default values of the bin (nBins parameter) is associated with the number of bins for human species, which depends on the length of the genome. nBins varies for different species, and will be set based on the input count matrix when running the simulation function. 



