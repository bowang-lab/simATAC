##PBMCs benchmarking###########################################################################################################

# Create cell by bin matrix from generated .h5 (snap) file
p.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/10x_PBMC_5k/10x_PBMC_5k.snap",
  sample="labels",
  num.cores=1
)
p.x.sp = addBmatToSnap(p.x.sp, bin.size=5000)

# Read true cell labels with cell barcodes list
metadata <- read.table('/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/10x_PBMC_5k/10x_PBMC_5k_metadata.tsv', header = TRUE)
metadata$barcode <- substring(metadata$barcode, first = 1, last = 16)
label <- sapply(p.x.sp@barcode,
                function(x) as.character(metadata[which(metadata$barcode == x),]$label))
p.x.sp@sample <- unlist(label)

# Define each cell group
index1 <- which(p.x.sp@sample == "1")
p.x.sp1 <- p.x.sp[index1 ,]
# p.x.sp1
index2  <- which(p.x.sp@sample == "2")
p.x.sp2 <- p.x.sp[index2,]
# p.x.sp2
index3 <- which(p.x.sp@sample == "3")
p.x.sp3 <- p.x.sp[index3,]
# p.x.sp3
index4 <- which(p.x.sp@sample == "4")
p.x.sp4 <- p.x.sp[index4,]
# p.x.sp4
index5 <- which(p.x.sp@sample == "5")
p.x.sp5 <- p.x.sp[index5,]
# p.x.sp5
index6 <- which(p.x.sp@sample == "6")
p.x.sp6 <- p.x.sp[index6,]
# p.x.sp6
index7 <- which(p.x.sp@sample == "7")
p.x.sp7 <- p.x.sp[index7,]
# p.x.sp7
index8 <- which(p.x.sp@sample == "8")
p.x.sp8 <- p.x.sp[index8,]
# p.x.sp8


# This function simulates each cell group from PBMCs separately, and saves all together in
# a .h5 file for performing evaluation.
# Inputs:
# version: An string to put all benchmarking plots and files for this simulation in a folder with
# that name.
# mean: The Gaussian noise mean to be used for simATAC simulation.
# sd: The Gaussian noise standard deviation to be used for simATAC simulation.
# Output: -
#
simulatePBMC <- function(version,  mean, sd){

  # Create a folder for saving simulated count matrix.
  dir.create(paste("Results/PBMCs/", version, sep = ""))
  id = paste("Results/PBMCs/", version, "/PBMCs", sep = "")
  gc()

  # Simulate each cell group with simATAC.
  sim2.1 <- simulate(p.x.sp1,  mean, sd, species = "human")
  gc()
  sim2.2 <- simulate(p.x.sp2,  mean, sd, species = "human")
  gc()
  sim2.3 <- simulate(p.x.sp3,  mean, sd, species = "human")
  gc()
  sim2.4 <- simulate(p.x.sp4,  mean, sd, species = "human")
  gc()
  sim2.5 <- simulate(p.x.sp5,  mean, sd, species = "human")
  gc()
  sim2.6 <- simulate(p.x.sp6,  mean, sd, species = "human")
  gc()
  sim2.7 <- simulate(p.x.sp7,  mean, sd, species = "human")
  gc()
  sim2.8 <- simulate(p.x.sp8,  mean, sd, species = "human")
  gc()

  # Combine simulated cell groups together to save for further analysis (for performing cell type clustering analysis).
  data <- rbind(t(assays(sim2.1)$counts), t(assays(sim2.2)$counts), t(assays(sim2.3)$counts), t(assays(sim2.4)$counts),
                t(assays(sim2.5)$counts), t(assays(sim2.6)$counts), t(assays(sim2.7)$counts), t(assays(sim2.8)$counts))
  data <- as(data, "dgCMatrix")
  label <- unlist(c(p.x.sp1@sample, p.x.sp2@sample, p.x.sp3@sample, p.x.sp4@sample,
                    p.x.sp5@sample, p.x.sp6@sample, p.x.sp7@sample, p.x.sp8@sample))
  gc()

  # Save simulated matrix with labels in a h5 file.

  if (file.exists(paste(id, "_sim_mat.h5", sep="")))
    #Delete file if it exists
    file.remove(paste(id, "_sim_mat.h5", sep=""))

  mat <- summary(data)
  h5createFile(paste(id, "_sim_mat.h5", sep=""))
  h5write(mat, paste(id, "_sim_mat.h5", sep=""), "sim")
  h5write(label, paste(id, "_sim_mat.h5", sep=""), "label")
  return()
}


# This function performs benchmarking on simATAC's simulated counts with a specific version
# by plotting each cell group's original and simulated parameters and comparing them, and extracting
# statistical similarity parameters.
# Input:
# version: An string to put all benchmarking plots and files for this simulation in a folder with
# that name.
# Output: -
#
benchmarkPBMC <- function(version){

  address <- "Results"
  name <- "PBMCs"

  # Load simulated cells for the input version parameter.
  count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                   ncol(p.x.sp@bmat),
                   nrow(p.x.sp@bmat))

  # Define cell groups for each label.
  cell.p.1 <- count[, which(colnames(count) == "1")]
  cell.p.2 <- count[, which(colnames(count) == "2")]
  cell.p.3 <- count[, which(colnames(count) == "3")]
  cell.p.4 <- count[, which(colnames(count) == "4")]
  cell.p.5 <- count[, which(colnames(count) == "5")]
  cell.p.6 <- count[, which(colnames(count) == "6")]
  cell.p.7 <- count[, which(colnames(count) == "7")]
  cell.p.8 <- count[, which(colnames(count) == "8")]

  # Perform benchmarking and plot simulated and real parameters for each cell group.
  plotFigures(t(p.x.sp1@bmat), cell.p.1, address, name, "Cell1", version)
  # plotFigures(t(p.x.sp2@bmat), cell.p.2, address, name, "Cell2", version)
  # plotFigures(t(p.x.sp3@bmat), cell.p.3, address, name, "Cell3", version)
  # plotFigures(t(p.x.sp4@bmat), cell.p.4, address, name, "Cell4", version)
  # plotFigures(t(p.x.sp5@bmat), cell.p.5, address, name, "Cell5", version)
  # plotFigures(t(p.x.sp6@bmat), cell.p.6, address, name, "Cell6", version)
  # plotFigures(t(p.x.sp7@bmat), cell.p.7, address, name, "Cell7", version)
  # plotFigures(t(p.x.sp8@bmat), cell.p.8, address, name, "Cell8", version)

  # # Perform benchmarking and plot simulated and real parameters for each cell group.
  # p.x.sp.sim2 <- p.x.sp
  # p.x.sp.sim2@bmat <- rbind(t(cell.p.1), t(cell.p.2), t(cell.p.3), t(cell.p.4),
  #                           t(cell.p.5), t(cell.p.6), t(cell.p.7), t(cell.p.8))
  # p.x.sp.sim2@barcode <- unlist(c(p.x.sp1@barcode, p.x.sp2@barcode, p.x.sp3@barcode, p.x.sp4@barcode,
  #                                 p.x.sp5@barcode, p.x.sp6@barcode, p.x.sp7@barcode, p.x.sp8@barcode))
  # p.x.sp.sim2@sample <- unlist(c(p.x.sp1@sample, p.x.sp2@sample, p.x.sp3@sample, p.x.sp4@sample,
  #                                p.x.sp5@sample, p.x.sp6@sample, p.x.sp7@sample, p.x.sp8@sample))
  # p.x.sp.sim2@file <- unlist(c(p.x.sp1@file, p.x.sp2@file, p.x.sp3@file, p.x.sp4@file,
  #                              p.x.sp5@file, p.x.sp6@file, p.x.sp7@file, p.x.sp8@file))
  # p.x.sp.sim2@metaData <- rbind(p.x.sp1@metaData, p.x.sp2@metaData, p.x.sp3@metaData, p.x.sp4@metaData,
  #                               p.x.sp5@metaData, p.x.sp6@metaData, p.x.sp7@metaData, p.x.sp8@metaData)
  # gc()
  #
  # # Cell type clustering with SnapATAC.
  # # Adjsut the number of reduced dimensions to 5
  # nmi <- SnapATACClustering(p.x.sp.sim2,
  #                           "human",
  #                           paste(address, "/", name, "/", version, "/", name, "_sim.txt", sep = ""))
  #
  # write(paste(name, version, as.character(nmi), sep = "    "),
  #       file=paste(address, "/", name, "/", name, "_clustering_metric.txt", sep = ""),
  #       append=TRUE)
  # print(paste(name, "clustering nmi:", nmi, sep = " "))
}


## Cluster real PBMCs cells with SnapATAC.
# Adjsut the number of reduced dimensions to 5
# name <- "PBMCs"
# address <- "Results"
# version <- "simATACV1_0"
# SnapATACClustering(p.x.sp,
#          "human",
#          paste(address, "/", name, "/", version, "/", name, "_real.txt", sep = ""))


# # initial seed: 60760
# simulatePBMC("simATACV1_0_0", 0, 0)
# simulatePBMC("simATACV1_-0.3_0.3", -0.3, 0.3)
# simulatePBMC("simATACV1_-0.4_0.4", -0.4, 0.4)
#
# benchmarkPBMC("simATACV1_0_0")
# benchmarkPBMC("simATACV1_-0.3_0.3")
# benchmarkPBMC("simATACV1_-0.4_0.4")


# # initial seed: 1000
# simulatePBMC("simATACV2_0_0", 0, 0)
# simulatePBMC("simATACV2_-0.3_0.3", -0.3, 0.3)
# simulatePBMC("simATACV2_-0.4_0.4", -0.4, 0.4)
#
# benchmarkPBMC("simATACV2_0_0")
# benchmarkPBMC("simATACV2_-0.3_0.3")
# benchmarkPBMC("simATACV2_-0.4_0.4")


# # initial seed: 10
# simulatePBMC("simATACV3_0_0", 0, 0)
# simulatePBMC("simATACV3_-0.3_0.3", -0.3, 0.3)
# simulatePBMC("simATACV3_-0.4_0.4", -0.4, 0.4)
#
# benchmarkPBMC("simATACV3_0_0")
# benchmarkPBMC("simATACV3_-0.3_0.3")
# benchmarkPBMC("simATACV3_-0.4_0.4")


