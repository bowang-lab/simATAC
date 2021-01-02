callPeak.new <- function(orig, sim, peak.num, bin.name, name){

  bin.mean.orig <- colSums(orig)/nrow(orig)
  bin.mean.sim <- colSums(sim)/nrow(sim)

  peak.index.orig <- order(bin.mean.orig, decreasing=TRUE)[1:peak.num]
  peak.index.sim <- order(bin.mean.sim, decreasing=TRUE)[1:peak.num]

  write(bin.name[peak.index.orig], paste(name, "_peaks_real.txt", sep = ""))
  write(bin.name[peak.index.sim], paste(name, "_peaks_simualted.txt", sep = ""))

  print(length(intersect(bin.name[peak.index.orig], bin.name[peak.index.sim]))/peak.num)
}



readSim <- function(file, row, col){

  sum <- h5read(file, "sim")
  label <- h5read(file, "label")
  count <- sparseMatrix(i=as.numeric(sum[,1]), j=as.numeric(sum[,2]), x = as.numeric(sum[,3]), dims = c(row, col))
  rownames(count) <- label
  return(count)
}


extract_peaks <- function (orig, name, version, bin.name, peak.num){

  sim <- readSim(paste("Results/", name, "/", version, "/", name, "_sim_mat.h5", sep = ""),
                 nrow(orig),
                 ncol(orig))

  print(dim(orig))
  print(dim(sim))

  callPeak.new(orig, sim, peak.num, bin.name, version)
}


run_all_version <- function(orig, name, version.name, bin.name, peak.num){
  extract_peaks(orig, name, paste(version.name, 'V1_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V2_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V3_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V4_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V5_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V6_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V7_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V8_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V9_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V10_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V11_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V12_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V13_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V14_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V15_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V16_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V17_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V18_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V19_0_0_1', sep=''), bin.name, peak.num)
  extract_peaks(orig, name, paste(version.name, 'V20_0_0_1', sep=''), bin.name, peak.num)
}




##Buenrostro2018

b.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Buenrostro_2018/Buenrostro2018_new.snap",
  sample="labels",
  do.par = TRUE,
  num.cores=1
)
b.x.sp = addBmatToSnap(b.x.sp, bin.size=5000)

# Read true cell labels with cell barcodes list.
metadata <- read.table('/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Buenrostro_2018/SnapATAC_metadata_Buenrostro_2018.tsv', header = TRUE)
label <- sapply(tolower(b.x.sp@barcode), function(x) as.character(metadata[which(tolower(metadata$file) == x),]$label))
b.x.sp@sample <- unlist(label)

# Remove unknown cell group.
index <- which(b.x.sp@sample != "UNK")
b.x.sp <- b.x.sp[index,]


##call peaks
name <- 'Buenrostro2018'
version.name <- 'Buenrostro'
peak.num = 30000
run_all_version(b.x.sp@bmat, name, version.name, bin.names.human, peak.num)




##Cusanovich2018
c.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/Cusanovich_2018_subset/Cusanovich2018_new.snap",
  sample="labels",
  num.cores=1
)
c.x.sp = addBmatToSnap(c.x.sp, bin.size=5000)

##call peaks
name <- 'Cusanovich2018'
version.name <- 'Cusanovich'
peak.num = 30000
run_all_version(c.x.sp@bmat, name, version.name, bin.names.mouse, peak.num)




##PBMCs
p.x.sp = createSnap(
  file="/home/zeinab/Documents/Zeinab/SIMATAC/data/benchmarking/10x_PBMC_5k/PBMCs_new.snap",
  sample="labels",
  num.cores=1
)
p.x.sp = addBmatToSnap(p.x.sp, bin.size=5000)

##call peaks
name <- 'PBMCs'
version.name <- 'PBMCs'
peak.num = 30000
run_all_version(p.x.sp@bmat, name, version.name, bin.names.human, peak.num)
