callPeak <- function(orig, sim, peak.num, name, title){
  bin.mean.orig <- colSums(orig)/nrow(orig)
  bin.mean.sim <- colSums(sim)/nrow(sim)

  names(bin.mean.orig) <- name
  names(bin.mean.sim) <- name

  indice.orig <- whichpart(bin.mean.orig, n=peak.num)
  indice.sim <- whichpart(bin.mean.sim, n=peak.num)

  write(names(indice.orig), paste(title, "_peaks_real.txt", sep = ""))
  write(names(indice.sim), paste(title, "_peaks_simualted.txt", sep = ""))

  com <- length(intersect(names(indice.orig), names(indice.sim)))

  print(com)
  print(com/5000)
}

whichpart <- function(x, n=30) {
  nx <- length(x)
  p <- nx-n
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}

file = "inst/extdata/Human_genome_coordinates.txt"
data <- read.table(file)
bin.names.human <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")

file = "inst/extdata/Mouse_genome_coordinates.txt"
data <- read.table(file)
bin.names.mouse <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")
########
address <- "/home/zeinab/Documents/Zeinab/simATAC/results"
name <- "Buenrostro"
version <- "simATACV1_0"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(b.x.sp.3@bmat),
                 nrow(b.x.sp.3@bmat))
callPeak(b.x.sp.3@bmat, t(count), 5000, bin.names.human, "Buenrostro")



address <- "/home/zeinab/Documents/Zeinab/simATAC/results"
name <- "PBMC"
version <- "simATACV1_0"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(p.x.sp.2@bmat),
                 nrow(p.x.sp.2@bmat))
callPeak(p.x.sp.2@bmat, t(count), 5000, bin.names.human, "PBMC")


address <- "/home/zeinab/Documents/Zeinab/simATAC/results"
name <- "Cusanovich"
version <- "simATACV1_0"
count <- readSim(paste(address, "/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 ncol(c.x.sp.1@bmat),
                 nrow(c.x.sp.1@bmat))
callPeak(c.x.sp.1@bmat, t(count), 5000, bin.names.mouse, "Cusanovich")



peak.data <- data.frame(first=c("Bulk peaks", "Bulk peaks", "Bulk peaks",
                                "Real peak bins", "Real peak bins", "Real peak bins",
                                "Simulated peak bins", "Simulated peak bins", "Simulated peak bins"),
                        second=c("Bulk peaks", "Real peak bins", "Simulated peak bins",
                                 "Bulk peaks", "Real peak bins", "Simulated peak bins",
                                 "Bulk peaks", "Real peak bins", "Simulated peak bins"),
                        per=c(1, 95326/97998, 95326/97998,
                              4812/4996, 1, 4827/4996,
                              4815/5000, 4827/5000, 1))
