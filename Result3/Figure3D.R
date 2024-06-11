#################
#/p300s/biosoft/app/python/python2023/envs/R4.2.0/bin/R
.libPaths("/xtdisk/jiapl_group/zhuxl/R/x86_64-pc-linux-gnu-library/4.2")
setwd("/p300s/jiapl_group/zhuxl/Microbiome/humann")
library(edgeR)
library(limma)

limm<-function(exp,control,case,sample,fdr,logfc){
  group1 <- rep('control',dim(control)[1])
  group2 <- rep('case',dim(case)[1])
  grouplist<-c(group1,group2)
  table(grouplist)
  design <- model.matrix(~0+factor(grouplist))
  colnames(design)=levels(factor(grouplist))
  names<-vector()
  samples<-as.vector(rbind(control,case)[,1])
  #for(j in 1:length(samples)){
  #a<-grep(samples[j],colnames(exp))
  #names<-c(names,a)}
  exprSet<-exp[,samples]
  dim(exprSet)
  rownames(design)=colnames(exprSet)

  fit <- lmFit(exprSet, design)
  cont.matrix <- makeContrasts("control-case", levels=design)
  fit2 <- contrasts.fit(fit,cont.matrix)
  fit2 <- eBayes(fit2)
  tempOutput<-topTable(fit2,coef=1, n=Inf)

  DEG_all<-data.frame(as.vector(rownames(tempOutput)),tempOutput$logFC,tempOutput$AveExpr,tempOutput$P.Value,tempOutput$adj.P.Val,stringsAsFactors = FALSE)
  colnames(DEG_all)<-c("Gene","logFC","AveExpr","pval","adjp")

  up_index<-which(DEG_all$pval <= fdr & DEG_all$logFC >= logfc)
  DEG_up<-DEG_all[up_index,]
  DEG_up$Deg<-"up"
  down_index<-which(DEG_all$pval <= fdr & DEG_all$logFC <= (-logfc))
  DEG_down<-DEG_all[down_index,]
  DEG_down$Deg<-"down"
  DEG_select<-rbind(DEG_up,DEG_down)
  dim(DEG_select)

  name1<-paste0(sample,".txt")
  name2<-paste0(sample,"_select.txt")
  write.table(DEG_all,name1,sep="\t",quote=FALSE,row.names =F)
  write.table(DEG_select,name2,sep="\t",quote=FALSE,row.names =F)
}



setwd("/p300s/jiapl_group/zhuxl/Microbiome/humann/Result_gene")

gene<-read.table("../Result_genefamilies.tsv.result",header=T,sep="\t",row.names=1)
colnames(gene)<-gsub("^X","",colnames(gene))
dim(gene) #57522 283
gene[1:6,1:6]
gene<-gene[rowSums(gene)!=0,colSums(gene)!=0]


group<-read.table("../SampleGroup",header=T,sep="\t")
length(intersect(colnames(gene),group[,1]))
need<-data.frame(sample=colnames(gene))
group<-merge(need,group,by="sample")
dim(group)
table(group$Group)



G1_gene<-data.frame(sample=group[group$Group=="G1",1])
non_G1_gene<-data.frame(sample=group[group$Group!="G1",1])
dim(G1_gene);dim(non_G1_gene)

G2_gene<-data.frame(sample=group[group$Group=="G2",1])
non_G2_gene<-data.frame(sample=group[group$Group!="G2",1])
dim(G2_gene);dim(non_G2_gene)

G3_gene<-data.frame(sample=group[group$Group=="G3",1])
non_G3_gene<-data.frame(sample=group[group$Group!="G3",1])
dim(G3_gene);dim(non_G3_gene)

G4_gene<-data.frame(sample=group[group$Group=="G4",1])
non_G4_gene<-data.frame(sample=group[group$Group!="G4",1])
dim(G4_gene);dim(non_G4_gene)

limm(gene,G1_gene,G2_gene,"exp_G1_VS_G2",0.05,1)
limm(gene,G1_gene,G3_gene,"exp_G1_VS_G3",0.05,1)
limm(gene,G1_gene,G4_gene,"exp_G1_VS_G4",0.05,1)
limm(gene,G2_gene,G3_gene,"exp_G2_VS_G3",0.05,1)
limm(gene,G2_gene,G4_gene,"exp_G2_VS_G4",0.05,1)
limm(gene,G3_gene,G4_gene,"exp_G3_VS_G4",0.05,1)

limm(gene,non_G1_gene,G1_gene,"exp_non_G1_VS_G1",0.05,1)
limm(gene,non_G2_gene,G2_gene,"exp_non_G2_VS_G2",0.05,1)
limm(gene,non_G3_gene,G3_gene,"exp_non_G3_VS_G3",0.05,1)
limm(gene,non_G4_gene,G4_gene,"exp_non_G4_VS_G4",0.05,1)

################
library(ggplot2)
library(ggrepel)
library(tidygraph)
library(dplyr)

file<-dir()
file<-file[grep("txt",file)]
file<-file[grep("LR",file)]


H<-data.frame()
for(j in 1:length(file)){
var<-file[j]
var1<-gsub("LR_","",var)
var1<-gsub(".txt","",var1)
data<-read.table(file[j],header=T,sep="\t")
head(data)
data$Group<-var1
H<-rbind(H,data)
}

data$FOLD<-rep("unchange",dim(data)[1])
data[which(data$logFC>1 & data$adjp < 0.05),]$FOLD<- "up"
data[which(data$logFC< -1 & data$adjp < 0.05),]$FOLD<- "down"
table(data$FOLD)




ta<-as.data.frame(table(data$FOLD))
up_num<-ta[which(ta[,1]=="up"),2]
down_num<-ta[which(ta[,1]=="down"),2]
unchange_num<-ta[which(ta[,1]=="unchange"),2]
data$FOLD<-as.factor(data$FOLD)
table(data$FOLD)

tapply(data$logFC,data$FOLD,summary)
data_posi<-data[which(data$logFC>1 & data$adjp < 0.05),]
data_posi = data_posi %>% group_by(FOLD) %>% top_n(n = 10, wt = logFC)
data_nege<-data[which(data$logFC< -1 & data$adjp < 0.05),]
data_nege = data_nege %>% group_by(FOLD) %>% top_n(n = -10, wt = logFC)
dim(data_posi);dim(data_nege)

p1<-ggplot()+
geom_point(data,mapping=aes(x=logFC, y=-log10(adjp),
color=FOLD))+
geom_vline(xintercept=c(-1,0,1),color="grey",linetype=2)+
geom_hline(yintercept=-log10(0.01),color="grey",linetype=2)+
xlab("log(FC)")+ylab("-log10(FDR)")+
theme_bw()+
ggtitle(var1)+
 theme(legend.position="top",
text=element_text(size=15),
plot.title=element_text(size=15,hjust=0.5),
#panel.border=element_blank(),
axis.line.x=element_line(color="black"))+
geom_point(aes(x=data_posi$logFC,y=-log10(data_posi$adjp)),color="pink")+
geom_point(aes(x=data_nege$logFC,y=-log10(data_nege$adjp)),color="lightblue")+
geom_point(aes(x=data_un$logFC,y=-log10(data_un$adjp)),color="grey")+
scale_color_manual(values=c("down"="lightblue","unchange"="grey","up"="pink"),
labels=c(paste0("down: ",down_num),paste0("unchange: ",unchange_num),paste0("up: ",up_num)))+
geom_label_repel(data=data_posi,aes(x=logFC, y=-log10(data_posi$adjp),label =Gene),
point.padding=unit(1, "lines"),max.overlaps=100,segment.colour = "grey",fill="pink", size=3,
fontface="bold",,color="black",arrow = arrow(length=unit(0.01,"npc")))+
geom_label_repel(data=data_nege,aes(x=logFC, y=-log10(data_nege$adjp),label =Gene),
point.padding=unit(1, "lines"),max.overlaps=100,segment.colour = "grey",fill="lightblue",size=3,
fontface="bold",color="black",arrow = arrow(length=unit(0.01,"npc")))


ggsave(paste0("Volcano_",var1,"_0.01_1.pdf"),p1,height=7,width=7)
ggsave(paste0("Volcano_",var1,"_0.01_1.png"),p1,height=7,width=7)
}

###########
library(devtools)
devtools::install_github("Proteomicslab57357/UniprotR")
library(UniprotR) 


file<-dir()
file<-file[grep("tsv",file)]

H<-data.frame()
for(j in 1:length(file)){
var<-file[j]
var1<-gsub("exp_","",var)
var1<-gsub("_label","",var1)
var1<-gsub(".tsv","",var1)
data<-read.table(file[j],header=T,sep="\t")
data$Group<-var1
H<-rbind(H,data)
}
table(H$Group,H$FOLD)
write.table(H,"UniRef90_Gene_DEG_Group.tsv",sep="\t",quote=F,row.names=F)
H1<-H[H$FOLD!="unchange",]
dim(H)
dim(H1)
write.table(H1,"UniRef90_Gene_DEG_Group_select.tsv",sep="\t",quote=F,row.names=F)


###############
library(openxlsx)
library(pheatmap)
library(ggplot2)


pathcoverage<-read.xlsx("../pathway.xlsx",sheet="coverage",rowNames=TRUE)
pathabundance<-read.xlsx("../pathway.xlsx",sheet="abundance",rowNames=TRUE)
colnames(pathcoverage)<-gsub("^X","",colnames(pathcoverage))
colnames(pathabundance)<-gsub("^X","",colnames(pathabundance))
dim(pathcoverage) #136 993
dim(pathabundance) #136 993
pathcoverage[1:6,1:6]
pathabundance[1:6,1:6]

pathcoverage1<-pathcoverage[rowSums(pathcoverage)!=0,colSums(pathcoverage)!=0]
pathabundance1<-pathabundance[rowSums(pathabundance)!=0,colSums(pathabundance)!=0]
dim(pathcoverage1) #136 93
dim(pathabundance1) #136 93

group<-read.table("../SampleGroup",header=T,sep="\t",row.names=1)
length(intersect(colnames(gene),group[,1]))

pdf("pathway_coverage.pdf",height=20,width=20)
png("pathway_coverage.png",height=2000,width=2000)
pheatmap(pathcoverage1,annotation_col=group,
  cluster_cols=T,cluster_rows=T,
  show_rownames=T,show_colnames=F,
  col=colorRampPalette(c("white","red"))(1000))
dev.off()


pdf("pathway_abundance.pdf",height=20,width=20)
png("pathway_abundance.png",height=2000,width=2000)
pheatmap(pathabundance1,annotation_col=group,
  cluster_cols=T,cluster_rows=T,
  show_rownames=T,show_colnames=F,
  col=colorRampPalette(c("white","red"))(1000))
dev.off()


var<-colnames(pathabundance1)
need<-data.frame(Sample=var,Group=group[var,])
head(need)
var<-as.data.frame(table(need$Group))

H<-data.frame()
for(i in 1:dim(var)[1]){
a1<-as.vector(var[i,1])
g1<-pathabundance1[,need[need$Group==a1,"Sample"]]
g1<- as.data.frame(unlist(apply(g1,1,function(x) table(x!=0))))
var2<-grep("TRUE",rownames(g1))
var2<-data.frame(path=rownames(g1)[var2],num=g1[var2,1])
var2$path<-gsub(".TRUE","",var2$path)
var2$Group<-a1
H<-rbind(H,var2)
}
library(reshape2)
H1<-dcast(H,path~Group,value.var="num")
H1[is.na(H1)]<-0
write.csv(H1,"pathabundance1_non0_stat_percentage.csv",row.names=F)

H1<-read.csv("pathabundance1_non0_stat_percentage.csv",row.names=1)
H1[,1]<-H1[,1]/var[1,2]
H1[,2]<-H1[,2]/var[2,2]
H1[,3]<-H1[,3]/var[3,2]
H1[,4]<-H1[,4]/var[4,2]

path1<-unique(c(rownames(H1[order(H1[,1],decreasing=T),])[1:15],
rownames(H1[order(H1[,2],decreasing=T),])[1:15],
rownames(H1[order(H1[,3],decreasing=T),])[1:15],
rownames(H1[order(H1[,4],decreasing=T),])[1:15]))

pathabundance2<-pathabundance1[path1,]
write.csv(pathabundance2,"pathabundance2_top10.csv")

pdf("pathway_abundance_top10.pdf",height=7,width=7)
#png("pathway_abundance_top10.png",height=700,width=700)
pheatmap(pathabundance2,annotation_col=group,
  cluster_cols=T,cluster_rows=T,
  show_rownames=T,show_colnames=F,
  cellheight=7,cellwidth=2,
  col=colorRampPalette(c("white","red"))(1000))
dev.off()

#####################

library(circlize)
mat1<-H1[path1,]
mat1<-as.matrix(mat1)
write.csv(mat1,"mat1.csv",quote=F)
#rownames(mat1)<-unlist(strsplit(path1,":"))[c(seq(2,length(path1)*2,2)-1)]


pdf("Ciros_Abundance_Pathway_Group_full.pdf",height=20,width=20)

pdf("mat1_example.pdf",height=7,width=7)
chordDiagram(
  mat1, #grid.col = grid.col, 
annotationTrack = "grid", 
  preAllocateTracks = list(
    track.height = 1#max(strwidth(unlist(dimnames(mat1))))
  )
)

circos.track(
  track.index = 1, panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter, CELL_META$ylim[1], 
      CELL_META$sector.index,  facing = "clockwise", 
      niceFacing = TRUE, adj = c(0, 0.5)
    )
  }, bg.border = NA
)

circos.clear()
dev.off()
