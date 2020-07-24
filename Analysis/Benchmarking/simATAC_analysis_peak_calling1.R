# This function writes the name of peak bins (of real and simATAC generated bin by cell matrices)
# with the highest averaged read counts to a file with the format of "chrX:start-end" for each bin
# in a line. There is a Python code (simATAC_analysis_peak_calling2.py) that reads the saved files
# and extracts the percentage of intersection between the list of peak regions from different sources
# (bulk peaks, real peak bins, simulated peak bins) for each pair of region list.                                                                                                                                                         peak bins) for each pair of region list.
# Inputs:
# orig.count: Original sparse bin by cell matrix.
# sim.count: Simulated sparse bin by cell matrix.
# peak.num: Number of peak bins to extract from original and simulated bin by cell matrices.
# bin.name: The bin names (different for species) with the format of "chrX:start-end" for each bin.
# name: The name of the benchmark dataset.
# Output: -
callPeak <- function(orig.count, sim.count, peak.num, bin.name, name){

  bin.mean.orig <- colSums(orig.count)/nrow(orig.count)
  bin.mean.sim <- colSums(sim.count)/nrow(sim.count)

  names(bin.mean.orig) <- bin.name
  names(bin.mean.sim) <- bin.name

  peak.index.orig <- order(bin.mean.orig, decreasing=TRUE)[1:peak.num]
  peak.index.sim <- order(bin.mean.sim, decreasing=TRUE)[1:peak.num]

  write(names(peak.index.orig), paste(name, "_peaks_real.txt", sep = ""))
  write(names(peak.index.sim), paste(name, "_peaks_simualted.txt", sep = ""))

}


# Read the bin names for human into a list
human.file = "Benchmarking_Data/Human_genome_coordinates.txt"
data <- read.table(human.file)
bin.names.human <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")

# Read the bin names for mouse into a list
mouse.file = "Benchmarking_Data/Mouse_genome_coordinates.txt"
data <- read.table(mouse.file)
bin.names.mouse <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")

# Perform peak calling on benchmarking datasets
address <- "../Benchmarking_Data"
version <- "simATACV1_0"
# Peak calling for Buenrostro2018 simualtion without noise

name <- "Buenrostro2018"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(b.x.sp@bmat),
                 nrow(b.x.sp@bmat))
callPeak(b.x.sp@bmat, t(count), 5000, bin.names.human, "Buenrostro")


# Peak calling for Cusanovich2018 simualtion without noise
name <- "Cusanovich2018"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(c.x.sp@bmat),
                 nrow(c.x.sp@bmat))
callPeak(c.x.sp@bmat, t(count), 5000, bin.names.mouse, "Cusanovich")


# Peak calling for PBMCs simualtion without noise
name <- "PBMCs"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(p.x.sp@bmat),
                 nrow(p.x.sp@bmat))
callPeak(p.x.sp@bmat, t(count), 5000, bin.names.human, "PBMC")
