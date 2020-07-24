# This function prints the average of clustering evaluation metrics, including normalized
# mutual information (nmi), adjusted mutual information (ami), and adjusted rand index (ari).
# Inputs:
# name: The name of the dataset ("Buenrostro2018", "Cusanovich2018", "PBMCs").
# version: A string indicating the version of the simulation.
# Output: -
analyze_clustering_metrics <- function(name, version = NULL){

  nmi <- read.table(paste("Results/", name, "/", name, "_clustering_metric.txt", sep = ""), header = FALSE, dec = ".")
  nmi <- nmi[which(nmi$V2 == version),]
  colnames(nmi) <- c("data", "version", "nmi")

  metrics <- read.table(paste("Results/", name, "/", name, "_clustering_ami_ari.txt", sep = ""), header = FALSE, dec = ".")
  colnames(metrics) <- c("data", "version", "ami", "ari")
  data <- merge(nmi, metrics, by = "version")

  print(paste("nmi:", mean(data$nmi)))
  print(paste("ami:", mean(data$ami)))
  print(paste("ari:", mean(data$ari)))
}

analyze_clustering_metrics("Buenrostro2018", c("simATACV1_0_0", "simATACV2_0_0", "simATACV3_0_0"))
analyze_clustering_metrics("Cusanovich2018", c("simATACV1_0_0", "simATACV2_0_0", "simATACV3_0_0"))
analyze_clustering_metrics("PBMCs", c("simATACV1_0_0", "simATACV2_0_0", "simATACV3_0_0"))

analyze_clustering_metrics("Buenrostro2018", c("simATACV1_-0.3_0.3", "simATACV2_-0.3_0.3", "simATACV3_-0.3_0.3"))
analyze_clustering_metrics("Cusanovich2018", c("simATACV1_-0.3_0.3", "simATACV2_-0.3_0.3", "simATACV3_-0.3_0.3"))
analyze_clustering_metrics("PBMCs", c("simATACV1_-0.3_0.3", "simATACV2_-0.3_0.3", "simATACV3_-0.3_0.3"))

analyze_clustering_metrics("Buenrostro2018", c("simATACV1_-0.4_0.4", "simATACV2_-0.4_0.4", "simATACV3_-0.4_0.4"))
analyze_clustering_metrics("Cusanovich2018", c("simATACV1_-0.4_0.4", "simATACV2_-0.4_0.4", "simATACV3_-0.4_0.4"))
analyze_clustering_metrics("PBMCs", c("simATACV1_-0.4_0.4", "simATACV2_-0.4_0.4", "simATACV3_-0.4_0.4"))
