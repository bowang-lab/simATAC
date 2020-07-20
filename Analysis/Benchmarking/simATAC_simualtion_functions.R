library(simATAC)
library(tictoc)


# This function gets a snap object (from SnapATAC package) containing a cell by bin matrix and
# simulates the same number of cells as given in the input, with specified input simulation parameters.
# Inputs:
# my.x.sp: An snap object containing cell by bin matrix in the bmat field.
# mean: The Gaussian noise mean to be used for simATAC simulation.
# sd: The Gaussian noise standard deviation to be used for simATAC simulation.
# species: The species of the input data, to be used for naming the bins in the simulation.
# Output:
# sim: A SingleCellExperiment object containing simulated parametes, returned by simATAC simulator.
#
simulate <- function(my.x.sp, mean, sd, species){
  tic("Estimation time:")
  object <- simATAC::simATACEstimate(t(my.x.sp@bmat))
  x <- toc()
  write(paste(as.character(nrow(my.x.sp@bmat)), x$toc-x$tic, sep = "     "), file = "Results/simATAC_estimation_time.txt", append = TRUE)

  object <- simATAC::setParameters(object,
                          nCells = nrow(my.x.sp@bmat),
                          species = species,
                          noise.mean = mean,
                          noise.sd = sd)

  tic("Simulation time:")
  sim <- simATAC::simATACSimulate(object)
  x <- toc()
  write(paste(as.character(nrow(my.x.sp@bmat)), x$toc-x$tic, sep = "     "), file = "Results/simATAC_simulation_time.txt", append = TRUE)

  return(sim)
}


# This function gets an stored simulated matrix and parameters as a .h5 file and returns a sparse bin by cell count matrix.
# Inputs:
# file: An .h5 file containing simulated counts in three columns (cell number, bin number, count value) in
# the "sim" header, with cell labels stored in the "label" header.
# row: Number of rows (bin)
# col: Number of columns (cell)
# Output:
# count: A sparse bin by cell matrix extracted from the input file, with cell labels as columns' name.
#
readSim <- function(file, row, col){

  sum <- h5read(file, "sim")
  label <- h5read(file, "label")
  count <- sparseMatrix(i=as.numeric(sum[,2]), j=as.numeric(sum[,1]), x = as.numeric(sum[,3]), dims = c(row, col))
  colnames(count) <- label
  return(count)
}


# This functions gets original and simulated sparse count matrices and plots their parameters and performs
# statistical comparison to evaluate the similarity of simulated parameters to the original one.
# Inputs:
# orig.count: Original sparse bin by cell matrix.
# sim.count: Simulated sparse bin by cell matrix.
# address: The address of the main folder that contains a folder for each of the benchmarking datasets.
# name: Name of the benchmarking dataset (Buenrostro2018, Cusanovich2018, PBMCs)
# type: Cell type of the input cell as string.
# version: The version of the simulation. For each version, there is a separate folder with version name in the
# benchmarking dataset folder.
# Output: -
#
plotFigures <- function(orig.count, sim.count, address, name, type, version){

  plotParams(orig.count, sim.count, address, name, type, version)
  EvaluateStats(orig.count, sim.count, address, name, type, version)
}

