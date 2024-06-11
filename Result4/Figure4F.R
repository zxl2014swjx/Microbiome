setwd("D:/Microbiome/RiskSCore")


library(LDlinkR)
library(gggrid)
library(EnsDb.Hsapiens.v75)
library(locuszoomr)
library(openxlsx)
library(qqman)


data<-read.table("Result_risk_score.risk_score.assoc.linear",header=T)
################


setwd("linear")

cand<-data[data$P<1e-5,]
dim(cand)
# cand<-data[data$P<5e-8,]
snp<-as.vector(cand$SNP)
snp<-snp[grep("rs",snp)]


for(i in 1:length(snp)){
snp1<-snp[i]
loc <- locus(data = data, chrom="CHR",pos="BP",p="P",
index_snp=snp1,, flank = 1e5,LD="r2",
ens_db = "EnsDb.Hsapiens.v75")
summary(loc)

#loc$xrange
#a<-as.data.frame(loc$TX);dim(a)
#b<-as.data.frame(loc$EX);dim(b)

pdf(paste0("Plot_",snp1,".pdf"),height=5,width=3)
locus_plot(loc,labels = "index",border=TRUE)
dev.off()
png(paste0("Plot_",snp1,".png"),height=400,width=250)
locus_plot(loc,labels = "index",border=TRUE)
dev.off()


}


#############

## full control

loc1 <- locus(data = data, chrom="CHR",pos="BP",p="P",
index_snp="rs3208970", flank = 1e5,LD="r2",
ens_db = "EnsDb.Hsapiens.v75")

loc2 <- locus(data = data, chrom="CHR",pos="BP",p="P",
index_snp="rs200509913", flank = 1e5,LD="r2",
ens_db = "EnsDb.Hsapiens.v75")

loc3 <-  locus(data = data, chrom="CHR",pos="BP",p="P",
index_snp="rs200724870",, flank = 1e5,LD="r2",
ens_db = "EnsDb.Hsapiens.v75")


pdf("res_sig_5e-8.pdf", width = 9, height = 4)
multi_layout(ncol = 3,
plots = {
locus_plot(loc1,axis.text=10,use_layout = FALSE,labels = "index",border=TRUE)
locus_plot(loc2,axis.text=10,use_layout = FALSE,labels = "index",border=TRUE)
locus_plot(loc3,axis.text=10,use_layout = FALSE,labels = "index",border=TRUE)
})
dev.off()






#############

genes <- c("STAT4", "IRF5", "UBE2L3")

# generate list of 'locus' class objects, one for each gene
loclist <- lapply(snp, locus,
                  data = data,
                  ens_db = "EnsDb.Hsapiens.v75",
                  LD = "r2")

## produce 3 locus plots, one on each page
pdf("myplot.pdf")
multi_layout(loclist)
dev.off()

## place 3 locus plots in a row on a single page
pdf("myplot.pdf")
multi_layout(loclist, ncol = 3, labels = "index")
dev.off()



