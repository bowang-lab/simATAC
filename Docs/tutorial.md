# simATAC Tutorial

Welcome to the simATAC simulator! 

simATAC is an R package developed for simulating single-cell ATAC sequencing (scATAC-seq) count matrices. simATAC is an easy to use simulation framework that given a group of cells having similar biological characteristics in the format of a bin by cell matrix as input, it generates a synthetic bin by cell matrix resembling the real samples. simATAC mainly performs two estimation and simulation steps to generate final counts. simATAC provides the offer to convert the simulated bin by cell matrix into the binary version, peak by cell, and any list of regions (as features) by cell matrices. There are separate functions for each step, and this tutorial gives an overview and introduction to simATAC’s functionality.

Assuming there is a scATAC-seq dataset with cells having similar biological characteristics (e.g. cell type), you can convert BAM files into a bin by cell (5 kbp window) matrix with any customized pipeline. We used Snaptools to generate .snap files which contain the bin by cell array. For running examples, we use the snap file of the [GSE99172](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99172) real scATAC-seq sample. We will skip the snap generation (See [here](https://github.com/r3fang/SnapATAC/wiki/FAQs#whatissnap) for how to generate a snap file). Instead, we will download the snap file, which is provided in the example folder. 

GSE99172 dataset includes 288 cells from chronic myelogenous leukemia cell type, and the generated bin by cell matrix (from snap file) contains cells having similar biological characteristics. You can also create the bin by cell matrix for multiple cell types, based on your original dataset. However, simATAC modelling is based on the assumption that input cells are biologically similar. You need to perform the simulation for each cell group separately. Given the GSE99172 dataset, we use simATAC to simulate 1000 or more scATAC-seq cells by learning read distribution across cells and bins.


## Table of Contents

- [Loading simATAC](#load)
- [simATACCount class](#simATACCount)
- [Set parameters](#set)
- [Get parameters](#get)
- [Estimation function](#estimation)
- [Simulation function](#simulation)


<a name="load"></a>**Loading simATAC**        
You need to load the simATAC R package to be able to use it.

```bash
> library(simATAC)
```

<a name="simATACCount"></a>**simATACCount class**        
We defined a simATACCount class, which is a class specifically designed for storing simATAC scATAC-seq simulation parameters. You can create a new simATACCount object by:

```bash
> object <- newsimATACCount()
```

You can print the default values for the simATACCount object parameters (except the non.zero.pro parameter, which is a large list of bins' non-zero cell proportions and will be updated during the simulation step):

```bash
> object
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
```

You can print the description for each parameter documented in the simATAC R package by running

```bash
> ?simATACCount
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

<a name="get"></a>**Get parameters**

To look at a particular parameter in the simATACCount object, we can run simATACget() function by giving the simATACCount object and the parameter's name as input:

```bash
> nBins <- simATACget(object, "nBins")
> nBins
[1] 642098
> species <- simATACget(object, "species")
> species
[1] "human"
> nCells <- simATACget(object, "nCells")
> nCells
[1] 500
```

To return a list of parameters, we can use getParameters() function by giving a list of parameters in the input:

```bash
> params <- getParameters(object, c("nBins", "species", "nCells"))
> params
$nBins
[1] 642098

$species
[1] "human"

$nCells
[1] 500
```

<a name="set"></a>**Set parameters**

We can manually adjust the parameters of the simATACCount object via setParameters() function by giving the object, and a list of parameters.

```bash
> object <- setParameters(object, nCells = 1000)
> object <- setParameters(object, noise.mean = -0.3, noise.sd = 0.3)
> object
An object of class "simATACCount"
Slot "nBins":
[1] 642098

Slot "nCells":
[1] 1000

Slot "seed":
[1] 16216

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
[1] -0.3

Slot "noise.sd":
[1] 0.3
```

<a name="estimation"></a>**Estimation function**

For each user-input, simATAC performs two core simulation steps: (i) estimating the model parameters based on the input bin by cell matrix, including the library sizes of the cells, the non-zero cell proportions of each bin and the average of reads per bin (bin mean); (ii) generating a bin by cell matrix that resembles the original input scATAC-seq data by sampling from Gaussian mixture and polynomial models with the estimated parameters. simATAC outputs a count matrix as a [`SingleCellExperiment`][SCE] object from SingleCellExperiment package, offering additional functions to convert it to other types of feature matrices.
 
simATAC allows us to estimate the parameters of the real count matrix by simATACEstimate() function. simATACEstimate() gets either a sparse matrix or a SingleCellExperiment object as input to perform the estimation. If the input matrix is in a SingleCellExperiment format, simATAC provides a getCountFromSCE() function that gets the SingleCellExperiment object as input and returns a sparse matrix as output. We can also directly pass the SingleCellExperiment object to the simATACEstimate() function. We first show an example of how to load the bin by cell matrix from the snap file.

simATAC provides a getCountFromh5() function, which is specifically implemented for reading the sparse bin by cell matrix created by Snaptools.

```bash
## return the cell by bin matrix from the snap file
> count <- getCountFromh5("GSE99172.snap")
## print dimentionality of the count matrix
> dim(count)
[1]    288 642098
```
The count object is a sparse matrix containing 288 cells, and 642098 bins with 5000 base-pair lengths. See the type of the count object by running

```bash
> typeof(count)
[1] "S4"
> class(count)
[1] "dgCMatrix"
attr(,"package")
[1] "Matrix"
```

The count object is a sparse cell by bin matrix (which is directly extracted from snap file), you need to convert it to a bin by cell matrix to be able to feed it into simATAC package.

```bash
> library(Matrix)
>
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

simATACEstimate function estimates the models' parameters, and if the verbose variable is TRUE (default value), it prints the progress of the estimation process. 
1. Library size parameters are estimated by fitting a Gaussian mixture model to the log-transformed of library size.
2. The proportion of cells having a non-zero count within each bin is calculated from the input count matrix (non-zero cell proportion).
3. The polynomial regression function parameters are estimated by fitting a quadratic function to the relationship between bin means and bin non-zero cell proportions.

All estimated model parameters are stored in the simATACCount object. 

The default value for the nBins parameter is associated with the number of bins for human species, which depends on the genome's length. nBins parameter varies for different species and will be adjusted based on the input count matrix when running the simATACSimulate() function. 

<a name="simulation"></a>**Simulation function**

Having the estimated parameters, simATAC generates the final counts by simATACSimulate() function. The number of cells to be simulated can be manually adjusted.

```bash
> sim <- simATACSimulate(object, nCells = 1000)
simATAC is:
...Updating parameters...
...Setting up SingleCellExperiment object...
...Simulating library size...
...Simulating non-zero cell proportion...
...Simulating bin mean...
...Generating final counts...
...Done...
```

```bash
> sim
class: SingleCellExperiment 
dim: 642098 1000 
metadata(1): Params
assays(1): counts
rownames(642098): chr1:1-5000 chr1:5001-10000 ...
  chrun_ki270392v1:1-5000 chrun_ki270394v1:1-5000
rowData names(2): Bin BinMean
colnames(1000): Cell1 Cell2 ... Cell999 Cell1000
colData names(2): Cell LibSize
reducedDimNames(0):
spikeNames(0):
altExpNames(0):
```

The simATACSimulate() function returns a sim object, a SingleCellExperiment object with 1000 cells in columns and 642098 (for human) bins stored in rows. You can access the simulated sparse count data via assays(), colData(), and rowData() functions from SummarizedExperiment R package. Row names of the SingleCellExperiment counts indicate the index of cells, and column names are associated with the positional information of the bins, including chromosome, starting position, and ending position.

```bash
> library(SingleCellExperiment)
> counts(sim)[1:5, 1:5]
5 x 5 sparse Matrix of class "dgCMatrix"
                 Cell1 Cell2 Cell3 Cell4 Cell5
chr1:1-5000          .     .     .     .     .
chr1:5001-10000      .     .     .     .     .
chr1:10001-15000     .     .     .     .     .
chr1:15001-20000     .     .     .     .     .
chr1:20001-25000     .     .     .     .     .
```

You can access the bin names and simulated bin means that are directly obtained from the polynomial function via rowData function from SummarizedExperiment package. Cells' names and library size values that are directly sampled from the Gaussian mixture model are also provided via colData function. Note that the BinMean and LibSize variables returned by colData and rowData are not from the final simulated counts, and as explained, they are the intermediary variables in the simulation process. 

```bash
> head(rowData(sim))
DataFrame with 6 rows and 2 columns
                              Bin             BinMean
                         <factor>           <numeric>
chr1:1-5000           chr1:1-5000 0.00282203544478947
chr1:5001-10000   chr1:5001-10000 0.00282203544478947
chr1:10001-15000 chr1:10001-15000 0.00407373696595591
chr1:15001-20000 chr1:15001-20000 0.00727215505384173
chr1:20001-25000 chr1:20001-25000 0.00282203544478947
chr1:25001-30000 chr1:25001-30000 0.00470551609382224
>
> head(colData(sim))
DataFrame with 6 rows and 2 columns
          Cell          LibSize
      <factor>        <numeric>
Cell1    Cell1 59456.8285774859
Cell2    Cell2 8431.56960395614
Cell3    Cell3 172995.376312487
Cell4    Cell4 68915.5581219072
Cell5    Cell5 22508.6642767295
Cell6    Cell6  43532.364826085
> 
```


[scater]: https://github.com/davismcc/scater
[SCE]: https://github.com/drisso/SingleCellExperiment