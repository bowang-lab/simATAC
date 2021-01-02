# This function prints the average of clustering evaluation metrics, including normalized
# mutual information (nmi), adjusted mutual information (ami), and adjusted rand index (ari).
# Inputs:
# name: The name of the dataset ("Buenrostro2018", "Cusanovich2018", "PBMCs").
# noise: A string indicating the noise level of the simulation.
# Output: -
analyze_clustering_metrics <- function(name, noise){

  data <- read.table(paste("Results/", name, "/", name, "_clustering_results.txt", sep = ""), header = FALSE, sep = ",")
  colnames(data) <- c("data", "version", "nmi", "ami", "ari")

  data <- data[grep(noise, data$version),]

  print(paste("nmi:", mean(data$nmi)))
  print(paste("ami:", mean(data$ami)))
  print(paste("ari:", mean(data$ari)))
}
