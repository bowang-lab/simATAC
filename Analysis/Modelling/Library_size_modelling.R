library(SnapATAC)
library(fitdistrplus)
library(mixtools)
library(diptest)
library(mclust)


# pdf("rplot_library_size.pdf")
# par(mfrow=c(3,3))
colname = c("dataset",
            "norm-chisqpvalue",
            "norm-kstest",
            "norm-ks",
            "diptest-pvalue",
            "n-ks-pvalue",
            "n-ks-D"
            )


results = list()
group <- vector()
# cellcount <- vector()
# libsize <- vector()


myFitDist <- function(lib.size, name){
  
  print(paste("**********",name,"**********"))
  
  ##mixtools package fitting result to multimodal distribution
  mixmdl = normalmixEM(lib.size, verb = FALSE)
  # plot(mixmdl, which=2)
  # lines(density(lib.size), lty=2, lwd=2)
  
  ##diptest package
  d <- dip(lib.size)
  dp.test <- dip.test(lib.size, simulate.p.value = FALSE, B = 2000)
  ##p-value computed via Monte Carlo simulation of a uniform distribution
  # dp.test <- dip.test(lib.size, simulate.p.value = TRUE, B = 2000)
  
  ##mclust package fitting result to unimodal distribution
  mod1 = densityMclust(lib.size, G = 1)
  
  ##mclust package fitting result to bimodal distribution
  mod2 <- densityMclust(lib.size, G = 2)
  # plot(mod2, what = "density", data = lib.size)
  # plot(mod2, what = "diagnostic", type = "cdf")
  # plot(mod2, what = "diagnostic", type = "qq")
  
  # library(lmtest)
  # lr <- lrtest (mod1, mod2)
  # lr.pval <- lr$`Pr(>Chisq)`[2]
  
  ##fitting normal distibution
  fn <- fitdist(lib.size, "norm")
  
  #Four Goodness-of-fit plots for normal distribution
  plot.legend <- c("normal")
  # denscomp(fn, legendtext = plot.legend)
  # qqcomp(fn, legendtext = plot.legend)
  # cdfcomp(fn, legendtext = plot.legend)
  # ppcomp(fn, legendtext = plot.legend)
  v <- gofstat(fn, fitnames = plot.legend)
  
  # two-sided Kolmogorov-Smirnov test - normal distribution
  ks.n <- ks.test(lib.size, "pnorm", m=unname(fn$estimate[1]), sd=unname(fn$estimate[2]))
  
  result <- list(name,
                 round(unname(v$chisqpvalue[1]), 6),
                 unname(v$kstest[1]),
                 round(unname(v$ks[1]), 5),
                 unname(dp.test$p.value),
                 round(unname(ks.n$p.value), 5),
                 round(unname(ks.n$statistic), 5)
                )
  
  return(result)
}



analyze_data <- function(x.sp, name){
  
  lib.size = log2(Matrix::rowSums(x.sp@bmat)+1)
  lib.size <- lib.size[lib.size > 0]
  lib.size <- lib.size[!is.na(lib.size)]
  
  # cellcount <<- rbind(cellcount, nrow(x.sp@bmat))
  # libsize <<- rbind(libsize, mean(log2(Matrix::rowSums(x.sp@bmat)+1)))
  
  group <<- rbind(group, name)
  
  ##original library size
  tmp <- myFitDist(lib.size, name)
  results <<- rbind(results, tmp)
}



## GSE99172
x.sp = createSnap(
  file="GSE99172.snap",
  sample="GSE99172",
  num.cores=1
)
print(paste("**********GSE99172**********"))
showBinSizes("GSE99172.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
analyze_data(x.sp, "GSE99172")



## GSE74310
x.sp = createSnap(
  file="GSE74310.snap",
  sample="labels",
  num.cores=1
)
print(paste("**********GSE74310**********"))
showBinSizes("GSE74310.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

mydata <- h5read("/home/znavidi/projects/def-wanglab/ATAC-seq-data/scATAC-seq/GSE74310/GSE74310_cell_by_bin_SnapATAC.h5", name = "Matrix")
x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(x.sp@sample == "acute myeloid leukemia-blast cell")
index2 <- which(x.sp@sample == "acute myeloid leukemia-leukmeia stem cell")
index3 <- which(x.sp@sample == "lymphoid primed multipotent progenitor")
index4 <- which(x.sp@sample == "Monocyte")

x.sp1 <- x.sp[index1,]
x.sp2 <- x.sp[index2,]
x.sp3 <- x.sp[index3,]
x.sp4 <- x.sp[index4,]

analyze_data(x.sp1, "GSE74310_acute myeloid leukemia-blast cell")
analyze_data(x.sp2, "GSE74310_acute myeloid leukemia-leukmeia stem cell")
analyze_data(x.sp3, "GSE74310_lymphoid primed multipotent progenitor")
analyze_data(x.sp4, "GSE74310_Monocyte")
analyze_data(x.sp, "GSE74310")



## GSE65360
x.sp = createSnap(
  file="GSE65360.snap",
  sample="labels",
  num.cores=1
)
print(paste("**********GSE65360**********"))
showBinSizes("GSE65360.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

mydata <- h5read("GSE65360_cell_by_bin_SnapATAC.h5", name = "Matrix")
x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(x.sp@sample == "promyelocytic leukemia cells")
index2 <- which(x.sp@sample == "human embryonic stem cell line")
index3 <- which(x.sp@sample == "fibroblasts")
index4 <- which(x.sp@sample == "lymphoblastoid cells")
index5 <- which(x.sp@sample == "chronic myeloid leukemia cells")
index6 <- which(x.sp@sample == "erythroleukemia cell line")

x.sp1 <- x.sp[index1,]
x.sp2 <- x.sp[index2,]
x.sp3 <- x.sp[index3,]
x.sp4 <- x.sp[index4,]
x.sp5 <- x.sp[index5,]
x.sp6 <- x.sp[index6,]

analyze_data(x.sp1, "GSE65360_promyelocytic leukemia cells")
analyze_data(x.sp2, "GSE65360_human embryonic stem cell line")
analyze_data(x.sp3, "GSE65360_fibroblasts")
analyze_data(x.sp4, "GSE65360_lymphoblastoid cells")
analyze_data(x.sp5, "GSE65360_chronic myeloid leukemia cells")
analyze_data(x.sp6, "GSE65360_erythroleukemia cell line")
analyze_data(x.sp, "GSE65360")



## SRR1947692
x.sp = createSnap(
  file="SRR1947692.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR1947692.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
print(paste("**********GSE68103/GM12878_HEK293T**********"))

mydata <- h5read("GM12878vsHEK_cell_by_bin_SnapATAC.h5", name = "Matrix")
bar.com <- intersect(x.sp@barcode, mydata$Label[,2])
index <- which(x.sp@barcode %in% bar.com)
new.x.sp <- x.sp[index,]
new.x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(new.x.sp@sample == "HEK293T")
index2 <- which(new.x.sp@sample == "GM12878")
index3 <- which(new.x.sp@sample == "Mixed")

x.sp1 <- new.x.sp[index1,]
x.sp2 <- new.x.sp[index2,]
x.sp3 <- new.x.sp[index3,]

analyze_data(x.sp1, "SRR1947692_HEK293T")
analyze_data(x.sp2, "SRR1947692_GM12878")
analyze_data(x.sp3, "SRR1947692_Mixed")
analyze_data(x.sp, "SRR1947692")



## SRR1947693
x.sp = createSnap(
  file="SRR1947693.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR1947693.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
print(paste("**********GSE68103/SRR1947693**********"))

mydata <- h5read("GM12878vsHL_cell_by_bin_SnapATAC.h5", name = "Matrix")
bar.com <- intersect(x.sp@barcode, mydata$Label[,2])
index <- which(x.sp@barcode %in% bar.com)
new.x.sp <- x.sp[index,]
new.x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(new.x.sp@sample == "HL60")
index2 <- which(new.x.sp@sample == "GM12878")
index3 <- which(new.x.sp@sample == "Mixed")

x.sp1 <- new.x.sp[index1,]
x.sp2 <- new.x.sp[index2,]
x.sp3 <- new.x.sp[index3,]

analyze_data(x.sp1, "SRR1947693_HL60")
analyze_data(x.sp2, "SRR1947693_GM12878")
analyze_data(x.sp3, "SRR1947693_Mixed")
analyze_data(x.sp, "SRR1947693")



## GSE112091
x.sp = createSnap(
  file="Breast_Tumor.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("Breast_Tumor.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
print(paste("**********GSE112091/GSE112245**********"))

mydata <- h5read("Breast_Tumor_SCALE_cell_by_bin_SnapATAC.h5", name = "Matrix")
tmp <- lapply(x.sp@barcode, tolower)
x.sp@barcode <- unlist(tmp)
bar.com <- intersect(x.sp@barcode, mydata$Label[,2])
index <- which(x.sp@barcode %in% bar.com)
new.x.sp <- x.sp[index,]
new.x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(new.x.sp@sample == "Epcam+")
index2 <- which(new.x.sp@sample == "CD45+")

x.sp1 <- new.x.sp[index1,]
x.sp2 <- new.x.sp[index2,]

analyze_data(x.sp1, "GSE112091_Epcam+")
analyze_data(x.sp2, "GSE112091_CD45+")
analyze_data(x.sp, "GSE112091/GSE112245")



## GSE100033
x.sp = createSnap(
  file="SRR6768122.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR6768122.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)
print(paste("**********GSE100033/GSM2668124**********"))

mydata <- h5read("SRR6768122_SCALE_cell_by_bin_SnapATAC.h5", name = "Matrix")
bar.com <- intersect(x.sp@barcode, mydata$Label[,2])
index <- which(x.sp@barcode %in% bar.com)
new.x.sp <- x.sp[index,]
new.x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(new.x.sp@sample == "EX3")
index2 <- which(new.x.sp@sample == "OC")
index3 <- which(new.x.sp@sample == "MG")
index4 <- which(new.x.sp@sample == "EX2")
index5 <- which(new.x.sp@sample == "IN2")
index6 <- which(new.x.sp@sample == "AC")
index7 <- which(new.x.sp@sample == "IN1")
index8 <- which(new.x.sp@sample == "EX1")

x.sp1 <- new.x.sp[index1,]
x.sp2 <- new.x.sp[index2,]
x.sp3 <- new.x.sp[index3,]
x.sp4 <- new.x.sp[index4,]
x.sp5 <- new.x.sp[index5,]
x.sp6 <- new.x.sp[index6,]
x.sp7 <- new.x.sp[index7,]
x.sp8 <- new.x.sp[index8,]

analyze_data(x.sp1, "GSE100033_GSM2668124_EX3")
analyze_data(x.sp2, "GSE100033_GSM2668124_OC")
analyze_data(x.sp3, "GSE100033_GSM2668124_MG")
analyze_data(x.sp4, "GSE100033_GSM2668124_EX2")
analyze_data(x.sp5, "GSE100033_GSM2668124_IN2")
analyze_data(x.sp6, "GSE100033_GSM2668124_AC")
analyze_data(x.sp7, "GSE100033_GSM2668124_IN1")
analyze_data(x.sp8, "GSE100033_GSM2668124_EX1")
analyze_data(x.sp, "GSE100033")



## GSE129785
sample.address <- list.files("/GSE129785/new_snap", pattern="*.snap$", all.files=FALSE, full.names=TRUE)
sample.name <- list.files("/GSE129785/new_snap", pattern="*.snap$", all.files=FALSE, full.names=FALSE)


for(i in 1:63){
  x.sp = createSnap(
    file=sample.address[i],
    sample=sample.name[i],
    num.cores=1
  )
  
  x.sp = addBmatToSnap(x.sp, bin.size=5000)
  print(paste(sample.name[i], " (", as.character(length(x.sp@barcode)), ")", sep = ""))
  analyze_data(x.sp, paste("GSE129785_", sample.name[i]))
}



colnames(results) <- colname
write.table(results, file = "Library_size_original.csv", append=F, sep='\t', row.names = TRUE, quote = FALSE, col.names = TRUE)
# dev.off()
