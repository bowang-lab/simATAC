#' Estimate and simulate simATAC simulation parameters
#'
#' Estimate parameters for the simATAC simulation from a real bin by cell
#' input matrix, and use them for simualting final counts.
#'
#' @param count Either a sparse bin by cell count matrix, or a SingleCellExperiment
#'        object containing the count matrix to estimate parameters from.
#' @param object A simATACCount object to store estimated parameters and
#'        the count matrix in it.
#' @param default Logical variable. Sets default parameters if TRUE.
#' @param verbose Logical variable. Prints the simulation progress if TRUE.
#' @param ... Any additional parameter settings to override what is provided in
#'        \code{simATACCount} object.
#'
#' @return SingleCellExperiment object containing the estimated counts and parameters.
#'
#'
#' @example
#' count <- getCountFromh5("GSE99172.snap")
#'
#' # simple simulation
#' sim <- simATACGenerate(coutn=t(count))
#'
#'\dontrun{
#' # set nCells parameter
#' sim <- simATACGenerate(count=count, nCells=500)
#'
#' # simulation with default parameters
#' sim <- simATACGenerate(default = TRUE, nCells=500)
#' }
#'
#' @seealso
#' \code{\link{simATACEstimate}}, \code{\link{setParameters}}
#' \code{\link{simATACSimulate}}
#' @export
#'
simATACGenerate <- function(count = NULL, object = newsimATACCount(), default = TRUE, verbose = TRUE, ...) {


  if(is.null(count)){
    default <- TRUE
  }
  if (default == FALSE){
    object <- simATACEstimate(count, object, verbose)
  }

  object <- setParameters(object, ...)
  sim <- simATACSimulate(object)

  return(sim)
}


#' Return a sparse count matrix from a SingleCellExperiment object. If count matrix is missing
#' a warning is printed and the first assay is returned.
#'
#' @param sce SingleCellExperiment input object containing counts.
#'
#' @return A sparse matrix containing counts.
#'
#' @example
#' \dontrun{
#' # gets cell by bin sparse matrix
#' count <- getCountFromh5("GSE99172.snap")
#'
#' # create SingleCellExperiment object with bin by cell matrix
#' sce <- SingleCellExperiment(assays = list(counts = t(count)))
#'
#' object <- simATACEstimate(sce)
#' sim <- simATACSimulate(object)
#' }
#' @export
#'
getCountFromSCE <- function(sce) {

  checkmate::assertClass(sce, "SingleCellExperiment")

  if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
    count <- SingleCellExperiment::counts(sce)
  } else {
    warning("counts assay is missing, using the first assay instead")
    count <- SummarizedExperiment::assay(sce)
  }

  return(count)
}


#' Simulate a bin's mean from its non-zero cell proportion from estimated 
#' second degree polynomial coefficients.
#'
#' @param x non-zero cell proportion of a bin.
#' @param c0 coefficient of x power 0.
#' @param c1 coefficient of x power 1.
#' @param c2 coefficient of x power 2.
#'
#' @return A dependent variable y, which is the simulated bin mean.
#'
simBinMeans <- function(x, c0, c1, c2){

  y=c0+c1*x+c2*x^2
  y = max(0, y)

  return(y)
}


#' Extract count matrix from h5 file
#'
#' @param file Input .h5 (or snap) file to read cell by bin matrix from it.
#'
#' @return A sparse cell by bin count matrix from .h5 file.
#'
#' @example
#' \dontrun{
#' count <- getCountFromh5("GSE99172.snap")
#' }
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom rhdf5 h5read
#' @export
#'
getCountFromh5 <- function(file){

  count.x <- h5read(file, "/AM/5000/idx")
  count.y <- h5read(file, "/AM/5000/idy")
  count.value <- h5read(file, "/AM/5000/count")
  count.barcodeLen <- h5read(file, "/FM/barcodeLen")
  count.binChrom <- h5read(file, "/AM/5000/binChrom")
  count.binStart <- h5read(file, "/AM/5000/binStart")

  nCells <- length(count.barcodeLen)
  nBins <- length(count.binChrom)

  bin.names <- sapply(seq(1, nBins, 1), function(x) paste(count.binChrom[x], count.binStart[x], sep = "_"))
  count <- Matrix::sparseMatrix(count.x, count.y, x = count.value, dims = c(nCells, nBins))

  return(count)
}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the binary version.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return Sparse matrix containing binary version of the input simulated bin by cell matrix.
#'
#' @example
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetBinary(sim)
#'
#' @export
#'
simATACgetBinary <- function(sim){

  count.bin <- BiocGenerics::counts(sim)
  count.bin[count.bin > 0] <- 1
  count.bin <- as(count.bin, "dgCMatrix")

  return(count.bin)

}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the
#' peak by cell matrix by extracting bins having the highest means.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param peak.num Number of peak bins to extract from the input bin by cell matrix.
#'
#' @return Sparse matrix containing peak.num bins with the highest bin means.
#'
#' @example
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetCellByPeak(sim)
#'
#' @importFrom Matrix rowSums
#' @export
#'
simATACgetCellByPeak <- function(sim, peak.num = 5000){

  bin.mean <- rowSums(BiocGenerics::counts(sim))/ ncol(BiocGenerics::counts(sim))
  peak.index <- order(bin.mean, decreasing=TRUE)[1:peak.num]
  count.peak <- as(BiocGenerics::counts(sim)[peak.index,], "dgCMatrix")

  return(count.peak)
}


#' Get the name of the input region from bed file in the format of chr:start-end as a string. Start
#' and end are coordinates of the beginning and end of a specific bin with 5000 length base pairs.
#'
#' @param region Input region in the format of [chr, start, end].
#' @param bin.name List of bin names of the raw bin by cell matrix.
#'
#' @return List of bin names that has intersection with the input region.
#'
getBin <- function(region, bin.name){

  int.div.start <- as.numeric(region[2])%/%5000
  int.div.end <- as.numeric(region[3])%/%5000
  remainder.end <- as.numeric(region[3])%%5000

  start <- int.div.start
  end <- if (remainder.end == 0) int.div.end-1 else int.div.end

  bed.name <- sapply(start:end, function(x)
    paste(region[1], ":", as.character(x*5000+1), "-", as.character((x+1)*5000), sep = ""))

  bed.index <- which(bin.name == bed.name)
  bed.index <- paste(as.character(bed.index), collapse = ':')

  return(bed.index)
}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the region by cell
#' matrix, with regions defined in the input BED file.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param file.bed BED file containing chromosome, start, and end positions of regions.
#'
#' @return Sparse matrix containing bins having intersection with the regions in the BED file.
#'
#' @example
#' \dontrun{
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetCellByRegion(sim, file.bed = "file.bed")
#' }
#'
#' @export
#'
simATACgetCellByRegion <- function(sim, file.bed){

  bin.name <- rownames(BiocGenerics::counts(sim))

  bed <- as.data.frame(read.table(file.bed, header = FALSE, sep="", stringsAsFactors=FALSE, quote=""))
  bed <- bed[,1:3]
  colnames(bed) <- c("chr", "start", "end")

  bed.index <- apply(bed, 1, getBin, bin.name = bin.name)
  bed.index <- paste(as.character(bed.index), collapse = ':')

  bed.index <- strsplit(bed.index,":")
  bed.index <- unique(unlist(bed.index))

  count <- BiocGenerics::counts(sim)[as.numeric(bed.index), ]

  return(count)
}


#' Compare the main parameters of simATAC's simulated bin by cell matrix to real
#' matrix, including library size, bin mean, non-zero cell proportion, and the
#' relation between bin means and non-zero cell proportions.
#'
#' @param sim Simulated bin by cell sparse matrix.
#' @param real Real bin by cell sparse matrix.
#' @param address Folder address to save plots.
#' @param name Name of the data.
#'
#' @return
#'
#' @importFrom ggpubr ggboxplot 
#' @importFrom ggplot2 ggplot ggtitle geom_point theme theme_bw xlab ylab element_blank element_text ylim xlim aes
#' @importFrom grDevices png dev.off
#' @importFrom utils write.table
#' @importFrom stats cor
#' @export
#'
simATACCompare <- function(sim, real, address, name){

  # Library size
  sim.lib.size <- log2(colSums(sim)+1)
  real.lib.size <- log2(colSums(real)+1)

  sim.lib.size.df <- data.frame(group = "Simulated", value = sim.lib.size)
  real.lib.size.df <- data.frame(group = "Real", value = real.lib.size)

  df <- rbind(real.lib.size.df, sim.lib.size.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("#00AFBB", "#FC4E07"),
                 add = "jitter", shape = "group") + ggtitle(name) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  
  png(paste(address, "/Library_size_boxplot.png", sep = ""))
  print(p+ylab("log2(library size)"))
  dev.off()
  
  
  # Bin sparsity (proportion of zero entries in each bin)
  real.bin.sparsity <- rowSums(real == 0)/ncol(real)
  sim.bin.sparsity <- rowSums(sim == 0)/ncol(sim)
  
  real.bin.sparsity.df <- data.frame(group = "Real", value = real.bin.sparsity)
  sim.bin.sparsity.df <- data.frame(group = "Simulated", value =sim.bin.sparsity)
  
  df <- rbind(real.bin.sparsity.df, sim.bin.sparsity.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("red", "blue"),
                 add = "jitter", shape = "group") + ggtitle(name) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  
  png(paste(address, "/Bin_sparsity_boxplot.png", sep = ""))
  print(p+ylab("Bin sparsity")+ xlab(NULL))
  dev.off()
  
  
  # Cell sparsity (proportion of zero entries in each cell)
  real.cell.sparsity <- colSums(real == 0)/nrow(real)
  sim.cell.sparsity <- colSums(sim == 0)/nrow(sim)
  
  real.cell.sparsity.df <- data.frame(group = "Real", value = real.cell.sparsity)
  sim.cell.sparsity.df <- data.frame(group = "Simulated", value =sim.cell.sparsity)
  
  df <- rbind(real.cell.sparsity.df, sim.cell.sparsity.df)
  p <- ggboxplot(df, x = "group", y = "value",
                 color = "group", palette =c("#CC79A7", "#009E73"),
                 add = "jitter", shape = "group") + ggtitle(name) +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    ) + theme(legend.position = "none")
  
  png(paste(address, "/Cell_sparsity_boxplot.png", sep = ""))
  print(p+ylab("Cell sparsity")+ xlab(NULL))
  dev.off()
  
  
  # Relation between bin means and non-zero cell proportions.
  sim.bin.mean <- rowSums(sim)/ncol(sim)
  real.bin.mean <- rowSums(real)/ncol(real)

  sim.nzp <- rowSums(sim != 0)/ncol(sim)
  real.nzp <- rowSums(real != 0)/ncol(real)
  
  max <- max(c(real.bin.mean[which(real.nzp <= 0.8)], sim.bin.mean[which(real.nzp <= 0.8)]))
  
  # Plot original non-zero cell proportion and bin mean scatter plot (including bins with non-zero proportion <= 0.8 in real data)
  df <- data.frame(nzp = real.nzp[which(real.nzp<=0.8)], mean = real.bin.mean[which(real.nzp<=0.8)])
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(paste('Real', name, sep = " ")) +
    ylim(0, max) + xlim(0, 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )
  
  png(paste(address, "/Real_bin_mean_and_nonzero_proportion.png", sep = ""))
  print(p + xlab("Bin non-zero proportion") + ylab("Bin mean"))
  dev.off()
  
  
  # Plot simulated non-zero cell proportion and bin mean scatter plot (including bins with non-zero proportion <= 0.8 in real data)
  df <- data.frame(nzp = sim.nzp[which(sim.nzp<=0.8)], mean = sim.bin.mean[which(sim.nzp<=0.8)])
  p <- ggplot(df, aes(x = nzp, y = mean)) +
    geom_point() +
    ggtitle(paste('Simulated', name, sep = " ")) +
    ylim(0, max) + xlim(0, 0.8) +
    theme_bw() +
    theme(
      plot.title = element_text(size=24, face="bold.italic", hjust = 0.5),
      axis.title.x = element_text(size=24),
      axis.title.y = element_text(size=24),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24)
    )

  png(paste(address, "/Simulated_bin_mean_and_nonzero_proportion.png", sep = ""))
  print(p + xlab("Bin non-zero proportion") + ylab("Bin mean"))
  dev.off()
  
  
  if (dim(sim)[2] == dim(real)[2]){
    # Save Pearson correlation between bin means and non-zero cell proportion parameters.
    write.table(c(paste('Name', 'Parameter', 'Pearson_correlation', sep = "    ")),
                  file = paste(address, "/", name, "_correlation.txt", sep = ""),
                  append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)

    # Save the bin mean Pearson correlation between real and simulated datasets (including bins with non-zero proportion <= 0.8 in real data).
    write.table(c(paste(name, "BinMeanCor<=0.8", as.character(cor(real.bin.mean[which(real.nzp <= 0.8)], sim.bin.mean[which(real.nzp <= 0.8)], method = c("pearson"))), sep = "    ")),
                file = paste(address, "/", name, "_correlation.txt", sep = ""),
                append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE)

    # Save the non-zero cell proportion Pearson correlation between real and simulated datasets.
    write.table(c(paste(name, "NZeroPropCor", as.character(cor(real.nzp, sim.nzp, method = c("pearson"))), sep = "    ")),
                file = paste(address, "/", name, "_correlation.txt", sep = ""),
                append = TRUE, sep = " ", col.names = FALSE, quote = FALSE, row.names = FALSE) 
  }

}
