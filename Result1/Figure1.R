library(BiocManager)
library(phyloseq)
library(microbiome)
library(ape)
library(ggplot2)
library(compositions)
library(rbiom)
library(openxlsx)
library(knitr)
library(permute)
library(lattice)
library(vegan) 
library(ggpubr)
library(rlang)
library(cli)
library(forcats)
library(tidyverse) 
library(permute)
library(Maaslin2)
library(ggthemes)


setwd("D:/Microbiome")
mpa_table<-read.xlsx("CASPMI_Braken.xlsx",sheet="Species")

dim(mpa_table)
rownames(mpa_table)<-as.vector(mpa_table[,1])
mpa_table<-mpa_table[,-c(1)]
mpa_table[1:4,1:4]
mpa_table_bak[1:4,1:4]

mpa_table_bak<-mpa_table

mpa_table[mpa_table<=0.01]<-0
write.table(mpa_table,"Filter_mpa_table.xls",sep="\t",quote=F)

alpha_diversity <- microbiome::alpha(mpa_table, 
index = c("richness_observed","chao1", #Richness
"diversity_shannon","diversity_gini_simpson","dominance_gini","diversity_inverse_simpson","diversity_coverage",
"evenness_camargo","evenness_pielou","evenness_simpson","evenness_evar","evenness_bulla",
"dominance_dbp","dominance_dmn","dominance_absolute","dominance_relative","dominance_simpson","dominance_core_abundance",
"rarity_log_modulo_skewness","rarity_low_abundance","rarity_rare_abundance"))
head(alpha_diversity)
dim(alpha_diversity)



BCI<-t(mpa_table)
row<-rownames(BCI)
col<-colnames(BCI)
pinpu<-matrix(rep(0,length(row)*length(col)),nrow=length(row))
rownames(pinpu)<-row
colnames(pinpu)<-col
for(i in 1:length(row)){
for(j in 1:length(col)){
pinpu[i,j]=as.numeric(as.vector(BCI[i,j]))
}
}
BCI<-pinpu
write.csv(pinpu,"BCI_pinpu.csv")
BCI<-BCI*100
BCI[1:4,1:4]


shannon<-vegan::diversity(BCI,'shannon')
simpson<-vegan::diversity(BCI,'simpson')
invsimpson<-vegan::diversity(BCI,'inv')
unbias.simpson <- vegan::simpson.unb(BCI)
alpha <- vegan::fisher.alpha(BCI)
species <- vegan::specnumber(BCI) 
evenness <- species/log(species) 
data<-data.frame(shannon,simpson,invsimpson,
unbias.simpson, alpha,species,evenness)

if(names(table(rownames(data)==rownames(alpha_diversity)))){
al<-cbind(data,alpha_diversity)
}
dim(al)
head(al[,c(1:6)])
ID<-rownames(al)
al<-cbind(ID,al)
write.table(al,"Microbiome_Diversity.xls",
sep="\t",quote=F,row.names=F)

#######################
head(mpa_table[,c(1:6)])
#mat <- rbiom::beta.div(as.matrix(mpa_table), method="unifrac")

mat_manhattan <- as.matrix(rbiom::beta.div(as.matrix(mpa_table),
method="manhattan"))
head(mat_manhattan[,1:6])
mat_euclidean <- as.matrix(rbiom::beta.div(as.matrix(mpa_table),
method="euclidean"))
head(mat_euclidean[,1:6])
mat_braycurtis <- as.matrix(rbiom::beta.div(as.matrix(mpa_table),
method="bray-curtis"))
head(mat_braycurtis[,1:6])
mat_jaccard <- as.matrix(rbiom::beta.div(as.matrix(mpa_table),
method="jaccard"))
head(mat_jaccard[,1:6])

library(pheatmap)
par(mfrow=c(2,2))
p1<-pheatmap(mat_manhattan,
show_rownames=F,show_colnames=F,main="manhattan dist")
p2<-pheatmap(mat_euclidean,
show_rownames=F,show_colnames=F,main="euclidean dist")
p3<-pheatmap(mat_braycurtis,
show_rownames=F,show_colnames=F,main="bray-curtis dist")
p4<-pheatmap(mat_jaccard,
show_rownames=F,show_colnames=F,main="jaccard dist")

plt<-cowplot::plot_grid(p1$gtable, p2$gtable, p3$gtable, p4$gtable,
ncol= 4)
ggsave("Beta_Diversity.pdf",plt,height=5,width=20)
ggsave("Beta_Diversity.png",plt,height=5,width=20)

#####################

phc <- function(d) { plot(hclust(d))}
phc(dist(mat_jaccard))

dis<-list(mat_manhattan,mat_euclidean,mat_braycurtis,mat_jaccard)
names(dis)<-c("mat_manhattan","mat_euclidean","mat_braycurtis","mat_jaccard")
li<-list()

for(i in 1:length(dis)){
adist<-as.dist(dis[[i]])
dist<-dis[[i]]
pc_num <-c(1,2)
pc_x <- pc_num[1]
pc_y <- pc_num[2]
pcoa <- cmdscale(dist, k=4, eig=TRUE)
pc12 <- pcoa$points[,pc_num]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits = 2)
pc12 <- as.data.frame(pc12)
colnames(pc12) <- c("pc_x","pc_y")
pc12['sample'] <- rownames(pc12)
head(pc12)

li[[i]]<-ggplot(pc12, aes(x = pc_x, y = pc_y))+
geom_point(alpha=0.5,)+
theme_bw()+
geom_vline(xintercept=0,linetype=2,color="orange")+
geom_hline(yintercept=0,linetype=2,color="blue")+
ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
ggtitle(paste0(names(dis)[i]," ",round(sum(pc[c(1,2)]),2),"%"))+
theme(plot.title = element_text(size=15,hjust = 0.5))
}


plt<-gridExtra::arrangeGrob(grobs=li,ncol=4)
ggsave("PCoA.pdf",plt,height=5,width=20)
ggsave("PCoA.png",plt,height=5,width=20)