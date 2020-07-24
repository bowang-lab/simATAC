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
[1] -0.1182035

Slot "mean.coef1":
[1] 10.34774

Slot "mean.coef2":
[1] -79.09698

Slot "noise.mean":
[1] 0

Slot "noise.sd":
[1] 0
and print the description for its parameters:

```

You can print the description for each parameter by running

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



