# This function prints the average of correlation parameters between real and simulated
# bin means, and real and simulated non-zero cell proportions (for each cell type).
# Three sets of simulation were perfomed for different sets of Gaussian noise parameters,
# and all simulations for a specific set were started from a specific seed value (first
# simulation set seed (version pattern = "simATACV1_mean_sd"): 60760, second simulation
# set seed (version pattern = "simATACV2_mean_sd"): 1000, third simulation set seed (version
# pattern = "simATACV3_mean_sd"): 10).
# This function reports the average of correlation of each parameter over all versions for
# simulations without Gaussian noise (mean:0, sd:0).
# Input:
# name: The name of the dataset ("Buenrostro2018", "Cusanovich2018", "PBMCs").
# Output: -
analyze_correlation <- function(name){

  data <- read.table(paste("Results/", name, "/", name, "_correlation.txt", sep = ""), header = FALSE, dec = ".")
  # returns correlation parameters only for simulation runs with simulation version equal to
  # "simATACV1_0_0", or "simATACV2_0_0", "simATACV3_0_0".
  data <- data[which("simATACV1_0_0" == data$V2 |
                       "simATACV2_0_0" == data$V2 |
                       "simATACV3_0_0" == data$V2),]

  # get all existing cell type or cell lable in the dataset and report the average of
  # correlation for each cell type separately.
  type <- unique(data$V1)

  for (t in type){
    print(paste(t,
                round(mean(data[which(data$V1 == t & data$V3 == "BinMeanPearCor<0.8"),]$V4), 2),
                round(mean(data[which(data$V1 == t & data$V3 == "NZeroPropPearCor"),]$V4), 2)
               )
        )
  }

}

analyze_correlation("Buenrostro2018")
analyze_correlation("Cusanovich2018")
analyze_correlation("PBMCs")
