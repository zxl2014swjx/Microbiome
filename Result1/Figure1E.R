setwd("D:/Microbiome/Combine/4Groups")

library(openxlsx)
library(grid)
library(futile.logger)
library(VennDiagram)
library(UpSetR)
library(dplyr)
library(tidyr)
require(ggplot2)
require(plyr)
require(gridExtra) 
require(grid)


data<-read.table("../mpa_table.xls",header=T,sep="\t")
dim(data)
data[1:4,1:4]
colnames(data)<-gsub("X","",colnames(data))


phenotype<-read.table("../mat_manhattan_pcoa_Group4.xls",header=T,sep="\t")
dim(phenotype)
head(phenotype)

G1<-as.vector(phenotype[phenotype$Group=="G4",]$sample)
length(G1)
n1<-data[,G1]
dim(n1)
rownames(n1)<-as.vector(rownames(data))
n1[1:4,1:4]
n2<-n1[rowSums(n1)!=0,]
n2[1:4,1:4]

s1<-rownames(n2)
s2<-rownames(n2)
s3<-rownames(n2)
s4<-rownames(n2)


length(s1);length(s2);length(s3);length(s4)

venn.diagram(list(Group1=s1,
Group2=s2,Group3=s3,Group4=s4),
width = 1300,height = 1300,
fill=c("red","green","blue","purple"), 
alpha=c(0.5,0.5,0.5,0.5), 
filename="Group1234_veen.tiff")

input <- c(
  'Group1' = length(s1),
  'Group2' = length(s2),
  'Group3' = length(s3),
  'Group4' = length(s4),
  'Group1&Group2' = length(intersect(s1,s2)),
  'Group1&Group3' = length(intersect(s1,s3)),
  'Group1&Group4' = length(intersect(s1,s4)),
  'Group2&Group3' = length(intersect(s2,s3)),
  'Group2&Group4' = length(intersect(s2,s4)),
  'Group3&Group4' = length(intersect(s3,s4)),
  'Group1&Group2&Group3' = length(intersect(intersect(s1,s2),s3)),
  'Group1&Group2&Group4' = length(intersect(intersect(s1,s2),s4)),
  'Group2&Group3&Group4' = length(intersect(intersect(s2,s3),s4)),
  'Group1&Group2&Group3&Group4' = length(intersect(intersect(intersect(s1,s2),s3),s4))
)

data <- fromExpression(input)

pdf("upset.pdf",height=5,width=7)
upset(data, nsets = 9, 
            sets = c('Group4','Group3',
                     'Group2' ,'Group1'),
            keep.order = TRUE, 
            point.size = 2, 
            line.size = 1.3, 
            mainbar.y.label = "IntersectionSize", 
            sets.x.label = "",
            mb.ratio = c(0.60, 0.40),
            text.scale = c(2, 2,2, 1.5,2, 2))
dev.off()


