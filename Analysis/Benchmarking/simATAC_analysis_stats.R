library(ggplot2)
library(scales)
library(reshape2)


# This function prints the average of Median Absolute Deviation (MAD), Mean Absolute Error
# (MAE) and Root Mean Square Error (RMSE) of each of the three main simualtion parameters,
# including library size, bin mean, and non-zero cell proportion, over all cell types in a
# dataset.
# Inputs:
# name: The name of the dataset ("Buenrostro2018", "Cusanovich2018", "PBMCs").
# version: A string indicating the version of the simulation.
# Output: -
analyze_MAD_MAE_RMSE_stats <- function(name, version){

  data <- read.table(paste("Results/", name, "/", name, "_stats_MAD_MAE_RMSE.txt", sep = ""), header = FALSE, dec = ".")
  data <- data[which(data$V2 == version),]

  data.mean <- data[which(data$V3 == "bin_mean"),]
  data.lib <- data[which(data$V3 == "lib_size"),]
  data.nzp <- data[which(data$V3 == "non_zero_proportion"),]

  mat <- matrix(0, nrow = 3, ncol = 3)
  mat[1,] <- c(mean(data.nzp[,4]), mean(data.nzp[,5]), mean(data.nzp[,6]))
  mat[2,] <- c(mean(data.mean[,4]), mean(data.mean[,5]), mean(data.mean[,6]))
  mat[3,] <- c(mean(data.lib[,4]), mean(data.lib[,5]), mean(data.lib[,6]))

  colnames(mat) <- c("MAD", "MAE", "RMSE")
  rownames(mat) <- c("Non-zero proportion", "Bin mean", "Library size")
  print(mat)
}

analyze_MAD_MAE_RMSE_stats("Buenrostro2018", "simATACV1_0_0")
analyze_MAD_MAE_RMSE_stats("Cusanovich2018", "simATACV1_0_0")
analyze_MAD_MAE_RMSE_stats("PBMCs", "simATACV1_0_0")
