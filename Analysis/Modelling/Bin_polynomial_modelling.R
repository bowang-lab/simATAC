library(SnapATAC)
library(logspline)
library(GenomicRanges)
library(outliers)
library(stats)
library(ggplot2)


# pdf("rplot_bin_distribution.png")
# par(mfrow=c(3,3))
r.squared <- vector()


my_fit_multinomial <- function(x, y, my_xlab, my_ylab, sign, name){
  
  fit1  <- lm(y~x)
  fit2 <- lm(y~poly(x, 2, raw=TRUE))
  # fit3 <- lm(y~poly(x, 3, raw=TRUE))
  
  sum <- summary(fit2)
  r.squared <<- rbind(r.squared, sum$r.squared)
  
  #   f1 <- function(x,c) coef(fit1)[1] + coef(fit1)[2]*x
  #   f2 <- function(x,c) coef(fit2)[1] + coef(fit2)[2]*x + coef(fit2)[3]*x^2
  #   # f3 <- function(x,c) coef(fit3)[1] + coef(fit3)[2]*x + coef(fit3)[3]*x^2 + coef(fit3)[4]*x^3
  #         
  #   df <- data.frame(x=x, y=y)
  #   pl <- ggplot(df) + geom_point(aes(x=x, y=y), size=2, colour="#993399") +
  #         xlab("Non-dropout probability") + ylab("Bin mean") +
  # 		    stat_function(fun=f1, size=1, show.legend = TRUE) +
  # 			  stat_function(fun=f2, size=1, show.legend = TRUE) +
  # 			  # stat_function(fun=f3, size=1, show.legend = TRUE) +
  # 			  ggtitle(paste(name, as.character(round(sum$r.squared, 2)), sep = "-")) +
  # 				theme( plot.title = element_text(hjust = 0.5, size = 12), text = element_text(size = 12)) +
  # 				scale_colour_manual("Function", labels = c("Degree 2", values = c("red")))
  # 	print(pl)
  
}



analyze_data <- function(x.sp, label){
  
  print(label)
  print(x.sp)
  
  cbb <- x.sp@bmat
  cbb <- cbb[which(rowSums(cbb) != 0),]
  
  lib.size <- rowSums(cbb)
  lib.med <- median(lib.size)
  norm.cbb <- cbb / lib.size * lib.med
  
  bin.mean <- colSums(norm.cbb)/nrow(norm.cbb)
  bin.nonzero.prob <- colSums(norm.cbb != 0)/nrow(norm.cbb)
  
  x <- bin.nonzero.prob[which(bin.nonzero.prob <= 0.8)]
  y <- bin.mean[which(bin.nonzero.prob <= 0.8)]
  
  my_fit_multinomial(x, y, "non-zero cells proportion <= 0.8", "bin mean", 1, label)
  
}



## GSE99172
x.sp = createSnap(
  file="GSE99172.snap",
  sample="label",
  num.cores=1
)
showBinSizes("GSE99172.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

mydata <- h5read("GSE99172_cell_by_bin_SnapATAC.h5", name = "Matrix")
x.sp@sample <- as.character(mydata$Label[,4])

index1 <- which(x.sp@sample == "chronic myelogenous leukemia")
x.sp1 <- x.sp[index1,]

analyze_data(x.sp1, "GSE99172")



## GSE74310
x.sp = createSnap(
  file="GSE74310.snap",
  sample="label",
  num.cores=1
)
showBinSizes("GSE74310.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

mydata <- h5read("GSE74310_cell_by_bin_SnapATAC.h5", name = "Matrix")
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



## GSE65360
x.sp = createSnap(
  file="GSE65360.snap",
  sample="label",
  num.cores=1
)
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



## RR1947692
x.sp= createSnap(
  file="SRR1947692.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR1947692.snap")
x.sp= addBmatToSnap(x.sp, bin.size=5000)

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



## SRR1947693
x.sp = createSnap(
  file="SRR1947693.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR1947693.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

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



## GSE112091
x.sp = createSnap(
  file="Breast_Tumor.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("Breast_Tumor.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

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



## GSE100033
x.sp = createSnap(
  file="SRR6768122.snap",
  sample="lable",
  num.cores=1
)
showBinSizes("SRR6768122.snap")
x.sp = addBmatToSnap(x.sp, bin.size=5000)

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



## GSE129785
snap <- list.files("GSE129785/new_snap", pattern="*.snap$", all.files=FALSE, full.names=TRUE)
names <- list.files("GSE129785/new_snap", pattern="*.snap$", all.files=FALSE, full.names=FALSE)

for (i in 1:length(snap)) {
  x.sp = createSnap(
    file=snap[i],
    sample="label",
    num.cores=1
  )
  x.sp = addBmatToSnap(x.sp, bin.size=5000)
  analyze_data(x.sp, paste("GSE129785_", substr(names[i],1,nchar(names[i])-5), sep = ""))
}


# dev.off()

