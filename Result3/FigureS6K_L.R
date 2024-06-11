
library(KEGGREST)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(openxlsx)
library(UniprotR) 

setwd("D:/Microbiome/3.humann/DEG_label")

data<-read.table("exp_G2_VS_G3_label.tsv",header=T,sep="\t")
dim(data)
head(data)

Accessions<-as.character(data[data$FOLD=="up",1])

TaxaObj <- GetNamesTaxa(Accessions,directorypath="./G2_VS_G3_up") 



PlotChromosomeInfo(TaxaObj,directorypath="./G2_VS_G3_up")
PlotGenesNetwork(TaxaObj,directorypath="./G2_VS_G3_up")

GeneOntologyObj <- GetProteinGOInfo(Accessions,directorypath="./G2_VS_G3_up") 
PlotGOBiological(GeneOntologyObj,directorypath="./G2_VS_G3_up") 
Plot.GOMolecular(GeneOntologyObj,directorypath="./G2_VS_G3_up")
Plot.GOSubCellular(GeneOntologyObj,directorypath="./G2_VS_G3_up")

PlotGoInfo(GeneOntologyObj,directorypath="./G2_VS_G3_up")


Pathway.Enr(Accessions,directorypath="./G2_VS_G3_up")

PlotEnrichedGO(Accs = Accessions, Path = "./G2_VS_G3_up")
PlotEnrichedPathways(Accs = Accessions, Path = "./G2_VS_G3_up")

PathologyObj <- GetPathology_Biotech(Accessions,directorypath="./G2_VS_G3_up")
Diseases <- Get.diseases(PathologyObj,directorypath="./G2_VS_G3_up")

GetproteinNetwork(Accessions,directorypath="./G2_VS_G3_up") 






