library(ggplot2)
library(ggpubr)
library(Matrix)
library(GenomicRanges)
library(SnapATAC)


# This function gets original and simulated bin by cell sparse matrices as input and extracts main simulation
# parameters (non-zero cell proportion, bin mean, and library size) and plot them to compare real and
# simulated counts.
# Inputs:
# orig.count: Original sparse bin by cell matrix.
# sim.count: Simulated sparse bin by cell matrix.
# address: The address of the main folder that contains a folder for each of the benchmarking datasets.
# name: Name of the benchmarking dataset (Buenrostro2018, Cusanovich2018, PBMCs)
# type: Cell type of the input cell as string.
# version: The version of the simulation. For each version, there is a separate folder with version name in the
# benchmarking dataset folder.
# Output: -
plotParams <- function(orig.count, sim.count, address, name, type, version){

  plot.title <- paste(name, "(", type, ")", sep = "")

  orig.nzp <- Matrix::rowSums(orig.count != 0)/ncol(orig.count)
  orig.mean <- Matrix::rowSums(orig.count)/ncol(orig.count)
  orig.libsize <- log2(Matrix::colSums(orig.count)+1)
  orig.bin.sparsity <- Matrix::rowSums(orig.count == 0)/ncol(orig.count)
  orig.cell.sparsity <- Matrix::colSums(orig.count == 0)/nrow(orig.count)

  sim.nzp <- Matrix::rowSums(sim.count != 0)/ncol(sim.count)
  sim.mean <- Matrix::rowSums(sim.count)/ncol(sim.count)
  sim.libsize <- log2(Matrix::colSums(sim.count)+1)
  sim.bin.sparsity <- Matrix::rowSums(sim.count == 0)/ncol(sim.count)
  sim.cell.sparsity <- Matrix::colSums(sim.count == 0)/nrow(sim.count)


  max <- max(c(orig.mean, sim.mean))

  # Plot original non-zero cell proportion and bin mean scatter plot (including all bins)
  df <- data.frame(nzp = orig.nzp, mean = orig.mean)
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(plot.title) +
    xlab("Bin non-zero proportion") +
    ylab("Bin mean") + ylim(0, max) + xlim(0, 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24, face="bold"),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_OrigBinMean&NZeroProp.png", sep = ""))
  print(p)
  dev.off()


  # Plot simulated non-zero cell proportion and bin mean relation (including all bins)
  df <- data.frame(nzp = sim.nzp, mean = sim.mean)
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(plot.title) +
    xlab("Bin non-zero proportion") +
    ylab("Bin mean") + ylim(0, max) + xlim(0, 1) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24, face="bold"),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_SimBinMean&NZeroProp.png", sep = ""))
  print(p)
  dev.off()


  max <- max(c(orig.mean[which(orig.nzp < 0.8)], sim.mean[which(orig.nzp < 0.8)]))

  # Plot original non-zero cell proportion and bin mean scatter plot (including bins with non-zero proportion < 0.8 in real data)
  df <- data.frame(nzp = orig.nzp[which(orig.nzp<0.8)], mean = orig.mean[which(orig.nzp<0.8)])
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(plot.title) +
    xlab("Bin non-zero proportion") +
    ylab("Bin mean") + ylim(0, max) + xlim(0, 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24, face="bold"),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_OrigBinMean&NZeroProp<0.8.png", sep = ""))
  print(p)
  dev.off()


  # Plot simulated non-zero cell proportion and bin mean scatter plot (including bins with non-zero proportion < 0.8 in real data)
  df <- data.frame(nzp = sim.nzp[which(sim.nzp<0.8)], mean = sim.mean[which(sim.nzp<0.8)])
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(plot.title) +
    xlab("Bin non-zero proportion") +
    ylab("Bin mean") + ylim(0, max) + xlim(0, 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24, face="bold"),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_SimBinMean&NZeroProp<0.8.png", sep = ""))
  print(p)
  dev.off()


  # Save the bin mean Pearson correlation between real and simulated datasets (including bins with non-zero proportion < 0.8 in real data).
  write.table(c(paste(plot.title, version, "BinMeanPearCor<0.8", as.character(cor(orig.mean[which(orig.nzp < 0.8)], sim.mean[which(orig.nzp < 0.8)], method = c("pearson"))), sep = "    ")),
              file = paste(address, "/", name, "/", name, "_correlation.txt", sep = ""),
              append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)


  # Save the non-zero cell proportion Pearson correlation between real and simulated datasets.
  write.table(c(paste(plot.title, version, "NZeroPropPearCor", as.character(cor(orig.nzp, sim.nzp, method = c("pearson"))), sep = "    ")),
              file = paste(address, "/", name, "/", name, "_correlation.txt", sep = ""),
              append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)


  orig.libsize.df <- data.frame(group = "Real", value = orig.libsize)
  sim.libsize.df <- data.frame(group = "Simulated", value = sim.libsize)

  # Plot box plot for original and simulated log transformed library sizes.
  df <- rbind(orig.libsize.df, sim.libsize.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("#00AFBB", "#FC4E07"),
                 add = "jitter", shape = "group") +
    ggtitle(plot.title) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_LibSizeBoxP.png", sep = ""))
  print(p+ylab("log2(library size)"))
  dev.off()


  orig.bin.sparsity.df <- data.frame(group = "Real", value = orig.bin.sparsity)
  sim.bin.sparsity.df <- data.frame(group = "Simulated", value =sim.bin.sparsity)

  # Plot box plot for original and simulated bin sparsity (proportion of zero entries in each bin).
  df <- rbind(orig.bin.sparsity.df, sim.bin.sparsity.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("red", "blue"),
                 add = "jitter", shape = "group") + ggtitle(plot.title) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_BinSparsityBoxP.png", sep = ""))
  print(p+ylab("Bin sparsity")+ xlab(NULL))
  dev.off()


  sx <- sort(orig.bin.sparsity); sy <- sort(sim.bin.sparsity)
  lenx <- length(sx)
  leny <- length(sy)
  if (leny < lenx)sx <- approx(1L:lenx, sx, n = leny)$y
  if (leny > lenx)sy <- approx(1L:leny, sy, n = lenx)$y
  df <- data.frame(orig = sx, sim = sy)
  p <- ggplot(df, aes(x = orig, y = sim)) +
    geom_point(color = "darkblue") +
    ggtitle(plot.title) +
    xlab("Real bin sparsity") +
    ylab("Simulated bin sparsity") +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24, face="bold"),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_", version, "_BinSparsityQQ.png", sep = ""))
  print(p + geom_abline(intercept = 0, slope = 1, col = "black"))
  dev.off()


  orig.cell.sparsity.df <- data.frame(group = "Real", value = orig.cell.sparsity)
  sim.cell.sparsity.df <- data.frame(group = "Simulated", value =sim.cell.sparsity)

  # Plot box plot for original and simulated cell sparsity (proportion of zero entries in each cell).
  df <- rbind(orig.cell.sparsity.df, sim.cell.sparsity.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("#CC79A7", "#009E73"),
                 add = "jitter", shape = "group") + ggtitle(plot.title) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24, face="bold"),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  # Save the plot in the given address
  png(paste(address, "/", name, "/", version, "/", name, "_", type, "_",  version, "_CellSparsityBoxP.png", sep = ""))
  print(p+ylab("Cell sparsity")+ xlab(NULL))
  dev.off()
}


# This function calculates the Median Absolute Deviation (MAD) between input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The MAD between real and simulated parameters.
getMAD <- function(orig, sim){
  mad <- median(abs(orig-sim))
  return(mad)
}


# This function calculates the Mean Absolute Error (MAE) between input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The MAE between real and simulated parameters.
getMAE <- function(orig, sim){
  mae <- mean(abs(orig-sim))
  return(mae)
}


# This function calculates the Root Mean Square Error (RMSE) between input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The RMSE between real and simulated parameters.
getRMSE <- function(orig, sim){
  # rmse <- sqrt(mean(abs(orig-sim)^2))
  rmse <- Metrics::rmse(orig, sim)
  return(rmse)
}


# This function calculates the Median Absolute Deviation (MAD) between sorted input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The MAD between sorted real and sorted simulated parameters.
getMAD.sorted <- function(orig, sim){
  mad <- median(abs(sort(orig)-sort(sim)))
  return(mad)
}


# This function calculates the Mean Absolute Error (MAE) between sorted input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The MAE between sorted real and sorted simulated parameters.
getMAE.sorted <- function(orig, sim){
  mae <- mean(abs((sort(orig)-sort(sim))))
  return(mae)
}


# This function calculates the Root Mean Square Error (RMSE) between sorted input vectors.
# Inputs:
# orig: Real data parameter.
# sim: Simulated data parameter.
# Output:
# mad: The RMSE between sorted real and sorted simulated parameters.
getRMSE.sorted <- function(orig, sim){
  # rmse <- sqrt(mean(abs((sort(orig)-sort(sim))^2)))
  rmse <- Metrics::rmse(sort(orig), sort(sim))
  return(rmse)
}


# This function calculates the Median Absolute Deviation (MAD), Mean Absolute Error (MAE), and
# Root Mean Square Error (RMSE) values for real and simulated main parameters: library size,
# bin mean, non-zero cell proportion.
# Inputs:
# orig.count: Original sparse bin by cell matrix.
# sim.count: Simulated sparse bin by cell matrix.
# address: The address of the main folder that contains a folder for each of the benchmarking datasets.
# name: Name of the benchmarking dataset (Buenrostro2018, Cusanovich2018, PBMCs)
# type: Cell type of the input cell as string.
# version: The version of the simulation. For each version, there is a separate folder with version name in the
# benchmarking dataset folder.
# Output: -
EvaluateStats <- function(orig.count, sim.count, address, name, type, version){

  ### library size ###
  orig.libsize <- log2(Matrix::colSums(orig.count)+1)
  sim.libsize <- log2(Matrix::colSums(sim.count)+1)

  mad <- getMAD.sorted(orig.libsize, sim.libsize)
  mae <- getMAE.sorted(orig.libsize, sim.libsize)
  rmse <- getRMSE.sorted(orig.libsize, sim.libsize)
  write.table(c(paste(paste(name,"(",type,")", sep = ""), version, "lib_size", as.character(mad), as.character(mae), as.character(rmse), sep = "    ")),
              file = paste(address, "/", name, "/", name, "_stats_MAD_MAE_RMSE.txt", sep = ""),
              append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)
  print(c(paste(name,"(",type,")", sep = ""), version, "lib", as.character(mad), as.character(mae), as.character(rmse)))

  ### bin means ###
  orig.binmean <- Matrix::rowSums(orig.count)/ncol(orig.count)
  sim.binmean <- Matrix::rowSums(sim.count)/ncol(sim.count)

  mad <- getMAD(orig.binmean, sim.binmean)
  mae <- getMAE(orig.binmean, sim.binmean)
  rmse <- getRMSE(orig.binmean, sim.binmean)
  write.table(c(paste(paste(name,"(",type,")", sep = ""), version, "bin_mean", as.character(mad), as.character(mae), as.character(rmse), sep = "    ")),
              file = paste(address, "/", name, "/", name, "_stats_MAD_MAE_RMSE.txt", sep = ""),
              append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)
  print(c(paste(name,"(",type,")", sep = ""), version, "mean", as.character(mad), as.character(mae), as.character(rmse)))

  ### bin non-dropout probability ###
  orig.non.drop <- Matrix::rowSums(orig.count != 0)/ncol(orig.count)
  sim.non.drop <- Matrix::rowSums(sim.count != 0)/ncol(sim.count)

  mad <- getMAD(orig.non.drop, sim.non.drop)
  mae <- getMAE(orig.non.drop, sim.non.drop)
  rmse <- getRMSE(orig.non.drop, sim.non.drop)
  write.table(c(paste(paste(name,"(",type,")", sep = ""), version, "non_zero_proportion", as.character(mad), as.character(mae), as.character(rmse), sep = "    ")),
              file = paste(address, "/", name, "/", name, "_stats_MAD_MAE_RMSE.txt", sep = ""),
              append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)
  print(c(paste(name,"(",type,")", sep = ""), version, "mean", as.character(mad), as.character(mae), as.character(rmse)))

}


# This function performs SnapATAC filteration and graph based clustering from SnapATAC website:
# (https://github.com/r3fang/SnapATAC)
# Inputs:
# x.sp: snap file.
# plot.file: File address to save
# Output: -
SnapATACClustering <- function(x.sp, species, label.file){

  x.sp <- x.sp[which(Matrix::rowSums(x.sp@bmat) != 0),]

  # Matrix binarization
  x.sp = makeBinary(x.sp, mat="bmat")

  # Bin filtration
  # First, filter out any bins overlapping with the ENCODE blacklist to prevent from potential artifacts.
  if (species == "mouse"){
    # system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz")
    black_list = read.table("mm10.blacklist.bed.gz")
  }else if(species == "human"){
    # system("wget http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/hg38-human/hg38.blacklist.bed.gz")
    black_list = read.table("hg38.blacklist.bed.gz")

    # OR download from https://www.encodeproject.org/files/ENCFF356LFX/
    # black_list = read.table("ENCFF356LFX.bed.gz")
  }

  black_list.gr = GRanges(
    black_list[,1],
    IRanges(black_list[,2], black_list[,3])
  )
  idy = queryHits(findOverlaps(x.sp@feature, black_list.gr))
  if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]};

  # Second, remove unwanted chromosomes.
  chr.exclude = seqlevels(x.sp@feature)[grep("random|chrM", seqlevels(x.sp@feature))]
  idy = grep(paste(chr.exclude, collapse="|"), x.sp@feature)
  if(length(idy) > 0){x.sp = x.sp[,-idy, mat="bmat"]}

  # Third, remove the top 5% bins that overlap with invariant features such as promoters of the house keeping genes.
  bin.cov = log10(Matrix::colSums(x.sp@bmat)+1)
  bin.cutoff = quantile(bin.cov[bin.cov > 0], 0.95)
  idy = which(bin.cov <= bin.cutoff & bin.cov > 0)
  x.sp = x.sp[, idy, mat="bmat"]
  x.sp

  # Dimensionality Reduction
  x.sp = runDiffusionMaps(
    obj=x.sp,
    input.mat="bmat",
    num.eigs=30
  )

  # Graph-based clustering
  x.sp = runKNN(
    obj=x.sp,
    eigs.dims=1:5,
    k=15
  )
  x.sp=runCluster(
    obj=x.sp,
    tmp.folder=tempdir(),
    louvain.lib="R-igraph",
    seed.use=10
  )
  x.sp@metaData$cluster = x.sp@cluster

  # Calculate Normalized Mutual Information (NMI) for predicted labels, and saveit in a file.
  library(aricode)
  true <- x.sp@sample
  pred <- x.sp@cluster
  nmi = NMI(true, pred)

  # Save true and predicted labels in a file.
  sample.cluster <- cbind(x.sp@sample, x.sp@cluster)
  colnames(sample.cluster) <- c("ground_truth", "SnapATAC_clusters")
  write.table(sample.cluster, label.file, append = FALSE, sep = "\t", dec = ".",
              row.names = FALSE, col.names = TRUE, quote = FALSE)

  return(nmi)
}
