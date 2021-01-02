library(SnapATAC)
library(tictoc)


# Create cell by bin matrix from generated .h5 (snap) file.
c.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Cusanovich_2018_subset/Cusanovich_2018_subset.snap",
  sample="labels",
  num.cores=1
)
c.x.sp = addBmatToSnap(c.x.sp, bin.size=5000)

# Read true cell labels with cell barcodes list.
metadata <- read.table('/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Cusanovich_2018_subset/Cusanovich_2018_subset_metadata.tsv', header = TRUE)
label <- sapply(c.x.sp@barcode,
                function(x) as.character(metadata[which(metadata$barcode == x),]$label))
c.x.sp@sample <- unlist(label)


# Define each cell group.
index1 <- which(c.x.sp@sample == "Lung")
c.x.sp1 <- c.x.sp[index1,]
# c.x.sp1
index2  <- which(c.x.sp@sample == "Thymus")
c.x.sp2 <- c.x.sp[index2,]
# c.x.sp2
index3 <- which(c.x.sp@sample == "Heart")
c.x.sp3 <- c.x.sp[index3,]
# c.x.sp3
index4 <- which(c.x.sp@sample == "Spleen")
c.x.sp4 <- c.x.sp[index4,]
# c.x.sp4
index5 <- which(c.x.sp@sample == "PreFrontalCortex")
c.x.sp5 <- c.x.sp[index5,]
# c.x.sp5
index6 <- which(c.x.sp@sample == "LargeIntestine")
c.x.sp6 <- c.x.sp[index6,]
# c.x.sp6
index7 <- which(c.x.sp@sample == "BoneMarrow")
c.x.sp7 <- c.x.sp[index7,]
# c.x.sp7
index8 <- which(c.x.sp@sample == "Liver")
c.x.sp8 <- c.x.sp[index8,]
# c.x.sp8
index9 <- which(c.x.sp@sample == "Cerebellum")
c.x.sp9 <- c.x.sp[index9,]
# c.x.sp9
index10 <- which(c.x.sp@sample == "SmallIntestine")
c.x.sp10 <- c.x.sp[index10,]
# c.x.sp10
index11 <- which(c.x.sp@sample == "Kidney")
c.x.sp11 <- c.x.sp[index11,]
# c.x.sp11
index12 <- which(c.x.sp@sample == "WholeBrain")
c.x.sp12 <- c.x.sp[index12,]
# c.x.sp12
index13 <- which(c.x.sp@sample == "Testes")
c.x.sp13 <- c.x.sp[index13,]
# c.x.sp13


# This function simulates each cell group from Cusanovich2018 separately, and saves all together in
# a .h5 file for performing evaluation.
# Inputs:
# version: An string to put all benchmarking plots and files for this simulation in a folder with
# that name.
# mean: The Gaussian noise mean to be used for simATAC simulation.
# sd: The Gaussian noise standard deviation to be used for simATAC simulation.
# Output: -
#
simulateCusanovich <- function(version, mean, sd){

  # Create a folder for saving simulated count matrix.
  dir.create(paste("Results/Cusanovich2018/", version, sep = ""))
  id = paste("Results/Cusanovich2018/", version, "/Cusanovich2018", sep = "")
  gc()

  # Simulate each cell group with simATAC.
  sim1.1 <- simulate(c.x.sp1,  mean, sd, species="mouse")
  gc()
  sim1.2 <- simulate(c.x.sp2,  mean, sd, species="mouse")
  gc()
  sim1.3 <- simulate(c.x.sp3,  mean, sd, species="mouse")
  gc()
  sim1.4 <- simulate(c.x.sp4,  mean, sd, species="mouse")
  gc()
  sim1.5 <- simulate(c.x.sp5,  mean, sd, species="mouse")
  gc()
  sim1.6 <- simulate(c.x.sp6,  mean, sd, species="mouse")
  gc()
  sim1.7 <- simulate(c.x.sp7,  mean, sd, species="mouse")
  gc()
  sim1.8 <- simulate(c.x.sp8,  mean, sd, species="mouse")
  gc()
  sim1.9 <- simulate(c.x.sp9,  mean, sd, species="mouse")
  gc()
  sim1.10 <- simulate(c.x.sp10,  mean, sd, species="mouse")
  gc()
  sim1.11 <- simulate(c.x.sp11,  mean, sd, species="mouse")
  gc()
  sim1.12 <- simulate(c.x.sp12,  mean, sd, species="mouse")
  gc()
  sim1.13 <- simulate(c.x.sp13,  mean, sd, species="mouse")
  gc()

  # Combine simulated cell groups together to save for further analysis (for performing cell type clustering analysis).
  data <- rbind(t(assays(sim1.1)$counts), t(assays(sim1.2)$counts), t(assays(sim1.3)$counts), t(assays(sim1.4)$counts),
                t(assays(sim1.5)$counts), t(assays(sim1.6)$counts), t(assays(sim1.7)$counts), t(assays(sim1.8)$counts),
                t(assays(sim1.9)$counts), t(assays(sim1.10)$counts), t(assays(sim1.11)$counts), t(assays(sim1.12)$counts),
                t(assays(sim1.13)$counts))
  data <- as(data, "dgCMatrix")
  label <- unlist(c(c.x.sp1@sample, c.x.sp2@sample, c.x.sp3@sample, c.x.sp4@sample,
                    c.x.sp5@sample, c.x.sp6@sample, c.x.sp7@sample, c.x.sp8@sample,
                    c.x.sp9@sample, c.x.sp10@sample, c.x.sp11@sample, c.x.sp12@sample,
                    c.x.sp13@sample))
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
benchmarkCusanovich <- function(version){

  address <- "Results"
  name <- "Cusanovich2018"

  # Load simulated cells for the input version parameter.
  count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                   ncol(c.x.sp@bmat),
                   nrow(c.x.sp@bmat))

  # Define cell groups for each label.
  cell.c.1 <- count[, which(colnames(count) == "Lung")]
  cell.c.2 <- count[, which(colnames(count) == "Thymus")]
  cell.c.3 <- count[, which(colnames(count) == "Heart")]
  cell.c.4 <- count[, which(colnames(count) == "Spleen")]
  cell.c.5 <- count[, which(colnames(count) == "PreFrontalCortex")]
  cell.c.6 <- count[, which(colnames(count) == "LargeIntestine")]
  cell.c.7 <- count[, which(colnames(count) == "BoneMarrow")]
  cell.c.8 <- count[, which(colnames(count) == "Liver")]
  cell.c.9 <- count[, which(colnames(count) == "Cerebellum")]
  cell.c.10 <- count[, which(colnames(count) == "SmallIntestine")]
  cell.c.11 <- count[, which(colnames(count) == "Kidney")]
  cell.c.12 <- count[, which(colnames(count) == "WholeBrain")]
  cell.c.13 <- count[, which(colnames(count) == "Testes")]

  # Perform benchmarking and plot simulated and real parameters for each cell group.
  plotFigures(t(c.x.sp1@bmat), cell.c.1, address, name, "Lung", version)
  plotFigures(t(c.x.sp2@bmat), cell.c.2, address, name, "Thymus", version)
  plotFigures(t(c.x.sp3@bmat), cell.c.3, address, name, "Heart", version)
  plotFigures(t(c.x.sp4@bmat), cell.c.4, address, name, "Spleen", version)
  plotFigures(t(c.x.sp5@bmat), cell.c.5, address, name, "PreFrontalCortex", version)
  plotFigures(t(c.x.sp6@bmat), cell.c.6, address, name, "LargeIntestine", version)
  plotFigures(t(c.x.sp7@bmat), cell.c.7, address, name, "BoneMarrow", version)
  plotFigures(t(c.x.sp8@bmat), cell.c.8, address, name, "Liver", version)
  plotFigures(t(c.x.sp9@bmat), cell.c.9, address, name, "Cerebellum", version)
  plotFigures(t(c.x.sp10@bmat), cell.c.10, address, name, "SmallIntestine", version)
  plotFigures(t(c.x.sp11@bmat), cell.c.11, address, name, "Kidney", version)
  plotFigures(t(c.x.sp12@bmat), cell.c.12, address, name, "WholeBrain", version)
  plotFigures(t(c.x.sp13@bmat), cell.c.13, address, name, "Testes", version)

  # Merge all cell groups together for cell type clustering analysis (with SnapATAC).
  c.x.sp.sim <- c.x.sp
  c.x.sp.sim@bmat <- rbind(t(cell.c.1), t(cell.c.2), t(cell.c.3), t(cell.c.4),
                            t(cell.c.5), t(cell.c.6), t(cell.c.7), t(cell.c.8),
                            t(cell.c.9), t(cell.c.10), t(cell.c.11), t(cell.c.12),
                            t(cell.c.13))
  c.x.sp.sim@barcode <- unlist(c(c.x.sp1@barcode, c.x.sp2@barcode, c.x.sp3@barcode, c.x.sp4@barcode,
                                 c.x.sp5@barcode, c.x.sp6@barcode, c.x.sp7@barcode, c.x.sp8@barcode,
                                 c.x.sp9@barcode, c.x.sp10@barcode, c.x.sp11@barcode, c.x.sp12@barcode,
                                 c.x.sp13@barcode))
  c.x.sp.sim@sample <- unlist(c(c.x.sp1@sample, c.x.sp2@sample, c.x.sp3@sample, c.x.sp4@sample,
                                c.x.sp5@sample, c.x.sp6@sample, c.x.sp7@sample, c.x.sp8@sample,
                                c.x.sp9@sample, c.x.sp10@sample, c.x.sp11@sample, c.x.sp12@sample,
                                c.x.sp13@sample))
  c.x.sp.sim@file <- unlist(c(c.x.sp1@file, c.x.sp2@file, c.x.sp3@file, c.x.sp4@file,
                              c.x.sp5@file, c.x.sp6@file, c.x.sp7@file, c.x.sp8@file,
                              c.x.sp9@file, c.x.sp10@file, c.x.sp11@file, c.x.sp12@file,
                              c.x.sp13@file))
  c.x.sp.sim@metaData <- rbind(c.x.sp1@metaData, c.x.sp2@metaData, c.x.sp3@metaData, c.x.sp4@metaData,
                               c.x.sp5@metaData, c.x.sp6@metaData, c.x.sp7@metaData, c.x.sp8@metaData,
                               c.x.sp9@metaData, c.x.sp10@metaData, c.x.sp11@metaData, c.x.sp12@metaData,
                               c.x.sp13@metaData)
  gc()

  # Cell type clustering with SnapATAC.
  # Adjsut the number of reduced dimensions to 5
  nmi <- SnapATACClustering(c.x.sp.sim,
                            "mouse",
                            paste(address, "/", name, "/", version, "/", name, "_sim.txt", sep = ""))

  write(paste(name, version, as.character(nmi), sep = "    "),
        file=paste(address, "/", name, "/", name, "_clustering_metric.txt", sep = ""),
        append=TRUE)
  print(paste(name, "clustering nmi:", nmi, sep = " "))

}

