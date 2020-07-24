library(SnapATAC)

# Create cell by bin matrix from generated .h5 (snap) file.
b.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Buenrostro_2018/Buenrostro_2018.snap",
  sample="labels",
  do.par = TRUE,
  num.cores=1
)
b.x.sp = addBmatToSnap(b.x.sp, bin.size=5000)

# Read true cell labels with cell barcodes list.
metadata <- read.table('/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Buenrostro_2018/SnapATAC_metadata_Buenrostro_2018.tsv', header = TRUE)
label <- sapply(b.x.sp@barcode,
              function(x) as.character(metadata[which(metadata$barcode == x),]$label))
b.x.sp@sample <- unlist(label)


# Remove unknown cell group.
index <- which(b.x.sp@sample != "UNK")
b.x.sp <- b.x.sp[index,]
# b.x.sp

# Define each cell group.
index1 <- which(b.x.sp@sample == "GMP")
b.x.sp1 <- b.x.sp[index1 ,]
# b.x.sp1
index2  <- which(b.x.sp@sample == "CMP")
b.x.sp2 <- b.x.sp[index2,]
# b.x.sp2
index3 <- which(b.x.sp@sample == "pDC")
b.x.sp3 <- b.x.sp[index3,]
# b.x.sp3
index4 <- which(b.x.sp@sample == "HSC")
b.x.sp4 <- b.x.sp[index4,]
# b.x.sp4
index5 <- which(b.x.sp@sample == "LMPP")
b.x.sp5 <- b.x.sp[index5,]
# b.x.sp5
index6 <- which(b.x.sp@sample == "MPP")
b.x.sp6 <- b.x.sp[index6,]
# b.x.sp6
index7 <- which(b.x.sp@sample == "mono")
b.x.sp7 <- b.x.sp[index7,]
# b.x.sp7
index8 <- which(b.x.sp@sample == "MEP")
b.x.sp8 <- b.x.sp[index8,]
# b.x.sp8
index9 <- which(b.x.sp@sample == "CLP")
b.x.sp9 <- b.x.sp[index9,]
# b.x.sp9


# This function simulates each cell group from Buenrostro2018 separately, and saves all together in
# a .h5 file for performing evaluation.
# Inputs:
# version: An string to put all benchmarking plots and files for this simulation in a folder with
# that name.
# mean: The Gaussian noise mean to be used for simATAC simulation.
# sd: The Gaussian noise standard deviation to be used for simATAC simulation.
# Output: -
#
simulateBuenrostro <- function(version, mean, sd){

  # Create a folder for saving simulated count matrix.
  dir.create(paste("Results/Buenrostro2018/", version, sep = ""))
  id = paste("Results/Buenrostro2018/", version, "/Buenrostro2018", sep = "")
  gc()

  # Simulate each cell group with simATAC.
  sim3.1 <- simulate(b.x.sp1, mean, sd, species = "human")
  gc()
  sim3.2 <- simulate(b.x.sp2, mean, sd, species = "human")
  gc()
  sim3.3 <- simulate(b.x.sp3, mean, sd, species = "human")
  gc()
  sim3.4 <- simulate(b.x.sp4, mean, sd, species = "human")
  gc()
  sim3.5 <- simulate(b.x.sp5, mean, sd, species = "human")
  gc()
  sim3.6 <- simulate(b.x.sp6, mean, sd, species = "human")
  gc()
  sim3.7 <- simulate(b.x.sp7, mean, sd, species = "human")
  gc()
  sim3.8 <- simulate(b.x.sp8, mean, sd, species = "human")
  gc()
  sim3.9 <- simulate(b.x.sp9, mean, sd, species = "human")
  gc()

  # Combine simulated cell groups together to save for further analysis (for performing cell type clustering analysis).
  data <- rbind(t(assays(sim3.1)$counts), t(assays(sim3.2)$counts), t(assays(sim3.3)$counts), t(assays(sim3.4)$counts),
                            t(assays(sim3.5)$counts), t(assays(sim3.6)$counts), t(assays(sim3.7)$counts), t(assays(sim3.8)$counts),
                            t(assays(sim3.9)$counts))
  data <- as(data, "dgCMatrix")
  label <- unlist(c(b.x.sp1@sample, b.x.sp2@sample, b.x.sp3@sample, b.x.sp4@sample,
                    b.x.sp5@sample, b.x.sp6@sample, b.x.sp7@sample, b.x.sp8@sample,
                    b.x.sp9@sample))
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
benchmarkBuenrostro <- function(version){

  address <- "Results"
  name <- "Buenrostro2018"

  # Load simulated cells for the input version parameter.
  count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                   ncol(b.x.sp@bmat),
                   nrow(b.x.sp@bmat))

  # Define cell groups for each label.
  cell.b.1 <- count[, which(colnames(count) == "GMP")]
  cell.b.2 <- count[, which(colnames(count) == "CMP")]
  cell.b.3 <- count[, which(colnames(count) == "pDC")]
  cell.b.4 <- count[, which(colnames(count) == "HSC")]
  cell.b.5 <- count[, which(colnames(count) == "LMPP")]
  cell.b.6 <- count[, which(colnames(count) == "MPP")]
  cell.b.7 <- count[, which(colnames(count) == "mono")]
  cell.b.8 <- count[, which(colnames(count) == "MEP")]
  cell.b.9 <- count[, which(colnames(count) == "CLP")]

  # Perform benchmarking and plot simulated and real parameters for each cell group.
  plotFigures(t(b.x.sp1@bmat), cell.b.1, address, name, "GMP", version)
  plotFigures(t(b.x.sp2@bmat), cell.b.2, address, name, "CMP", version)
  plotFigures(t(b.x.sp3@bmat), cell.b.3, address, name, "pDC", version)
  plotFigures(t(b.x.sp4@bmat), cell.b.4, address, name, "HSC", version)
  plotFigures(t(b.x.sp5@bmat), cell.b.5, address, name, "LMPP", version)
  plotFigures(t(b.x.sp6@bmat), cell.b.6, address, name, "MPP", version)
  plotFigures(t(b.x.sp7@bmat), cell.b.7, address, name, "mono", version)
  plotFigures(t(b.x.sp8@bmat), cell.b.8, address, name, "MEP", version)
  plotFigures(t(b.x.sp9@bmat), cell.b.9, address, name, "CLP", version)

  # Perform benchmarking and plot simulated and real parameters for each cell group.
  b.x.sp.sim <- b.x.sp
  b.x.sp.sim@bmat <- rbind(t(cell.b.1), t(cell.b.2), t(cell.b.3), t(cell.b.4),
                            t(cell.b.5), t(cell.b.6), t(cell.b.7), t(cell.b.8),
                            t(cell.b.9))
  b.x.sp.sim@barcode <- unlist(c(b.x.sp1@barcode, b.x.sp2@barcode, b.x.sp3@barcode, b.x.sp4@barcode,
                                  b.x.sp5@barcode, b.x.sp6@barcode, b.x.sp7@barcode, b.x.sp8@barcode,
                                  b.x.sp9@barcode))
  b.x.sp.sim@sample <- unlist(c(b.x.sp1@sample, b.x.sp2@sample, b.x.sp3@sample, b.x.sp4@sample,
                                 b.x.sp5@sample, b.x.sp6@sample, b.x.sp7@sample, b.x.sp8@sample,
                                 b.x.sp9@sample))
  b.x.sp.sim@file <- unlist(c(b.x.sp1@file, b.x.sp2@file, b.x.sp3@file, b.x.sp4@file,
                               b.x.sp5@file, b.x.sp6@file, b.x.sp7@file, b.x.sp8@file,
                               b.x.sp9@file))
  b.x.sp.sim@metaData <- rbind(b.x.sp1@metaData, b.x.sp2@metaData, b.x.sp3@metaData, b.x.sp4@metaData,
                                b.x.sp5@metaData, b.x.sp6@metaData, b.x.sp7@metaData, b.x.sp8@metaData,
                                b.x.sp9@metaData)
  gc()

  # Cell type clustering with SnapATAC.
  # Adjsut the number of reduced dimensions to 10
  nmi <- SnapATACClustering(b.x.sp.sim,
           "human",
           paste(address, "/", name, "/", version, "/", name, "_sim.txt", sep = ""))

  write(paste(name, version, as.character(nmi), sep = "    "),
        file=paste(address, "/", name, "/", name, "_clustering_metric.txt", sep = ""),
        append=TRUE)
  print(paste(name, "clustering nmi:", nmi, sep = " "))

}


## Cluster real Buenrostro2018 cells with SnapATAC.
# Adjsut the number of reduced dimensions to 10.
# name <- "Buenrostro"
# address <- "Results"
# version <- "simATACV1_0"
# analysis(b.x.sp,
#          "human",
#          paste(address, "/", name, "/", version, "/", name, "_real.txt", sep = ""))


# # initial seed: 60760
# simulateBuenrostro("simATACV1_0_0", 0, 0)
# simulateBuenrostro("simATACV1_-0.3_0.3", -0.3, 0.3)
# simulateBuenrostro("simATACV1_-0.4_0.4", -0.4, 0.4)
#
# benchmarkBuenrostro("simATACV1_0_0")
# benchmarkBuenrostro("simATACV1_-0.3_0.3")
# benchmarkBuenrostro("simATACV1_-0.4_0.4")


# # initial seed: 1000
# simulateBuenrostro("simATACV2_0_0", 0, 0)
# simulateBuenrostro("simATACV2_-0.3_0.3", -0.3, 0.3)
# simulateBuenrostro("simATACV2_-0.4_0.4", -0.4, 0.4)
#
# benchmarkBuenrostro("simATACV2_0_0")
# benchmarkBuenrostro("simATACV2_-0.3_0.3")
# benchmarkBuenrostro("simATACV2_-0.4_0.4")


# # initial seed: 10
# simulateBuenrostro("simATACV3_0_0", 0, 0)
# simulateBuenrostro("simATACV3_-0.3_0.3", -0.3, 0.3)
# simulateBuenrostro("simATACV3_-0.4_0.4", -0.4, 0.4)
#
# benchmarkBuenrostro("simATACV3_0_0")
# benchmarkBuenrostro("simATACV3_-0.3_0.3")
# benchmarkBuenrostro("simATACV3_-0.4_0.4")
