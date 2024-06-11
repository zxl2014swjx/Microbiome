setwd("D:/Microbiome/3.Gene")
library(DESeq2)

input_data<-read.table("All_Gene_Score.xls",header=T,sep="\t",row.names=1)
colnames(input_data)<-gsub("X","",colnames(input_data))
input_data <- round(input_data,digits = 0) 
input_data <- as.matrix(input_data)
input_data[1:4,1:4]
sample<-data.frame(Inx=1:dim(input_data)[2],SampleID=colnames(input_data))
anno<-read.table("meta.txt",header=T,sep="\t")
need<-merge(sample,anno)
dim(need)
need<-need[order(need$Inx),]
table(need$SampleID==colnames(input_data))
table(need$Group)

dim(input_data)
dim(need)

head(need)

group<-c("G1","G2","G3","G4")
table(need$Group)

need$G1<-rep("Neg",dim(need)[1])
need[need$Group=="G1",]$G1<-"Pos"
table(need$G1)

need$G2<-rep("Neg",dim(need)[1])
need[need$Group=="G2",]$G2<-"Pos"
table(need$G2)

need$G3<-rep("Neg",dim(need)[1])
need[need$Group=="G3",]$G3<-"Pos"
table(need$G3)

need$G4<-rep("Neg",dim(need)[1])
need[need$Group=="G4",]$G4<-"Pos"
table(need$G4)


H<-data.frame()
for(i in 4:dim(need)[2]){

condition <-factor(need[,i])
coldata <- data.frame(row.names=colnames(input_data),condition)

dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata,design=~condition) 
dds <- dds[rowSums(counts(dds))>1,] 
dds <-DESeq(dds)
res <- results(dds) 
summary(res)
head(dds)
res1<-as.data.frame(res)
res1$Gene<-rownames(res1)
head(res1)
res1<-na.omit(res1)

res1$Compare<-paste0(colnames(need)[i],"_VS_Other")

H<-rbind(H,res1)
}

table(H$Compare)

write.csv(H,)

colnames(H)[2]<-"avg_log2FC"
colnames(H)[5]<-"p_val"
colnames(H)[6]<-"p_val_adj"
colnames(H)[7]<-"gene"
colnames(H)[8]<-"cluster"
H$cluster<-gsub("_VS_Other","",H$cluster)
write.csv(H,"DESeq2_diff_gene.csv",row.names=F)

############


li<-list()
for(i in 1:4){


var1<-unique(H$cluster)[i]
data<-H[H$cluster==var1,]
data$DEG<-rep("unchange",dim(data)[1])
data[which(data$avg_log2FC>1 & data$p_val_adj < 0.05),]$DEG<- "up"
data[which(data$avg_log2FC< -1 & data$p_val_adj < 0.05),]$DEG<- "down"
table(data$DEG)
ta<-as.data.frame(table(data$DEG))
up_num<-ta[which(ta[,1]=="up"),2]
down_num<-ta[which(ta[,1]=="down"),2]
unchange_num<-ta[which(ta[,1]=="unchange"),2]
data$DEG<-as.factor(data$DEG)
table(data$DEG)
data_posi<-data[which(data$avg_log2FC>1 & data$p_val_adj < 0.05),]
data_nege<-data[which(data$avg_log2FC< -1 & data$p_val_adj < 0.05),]
data_un<-data[-which((data$avg_log2FC< -1 | data$avg_log2FC>1)& data$p_val_adj < 0.05),]
dim(data_posi);dim(data_nege);dim(data_un)

li[[i]]<-ggplot()+
geom_point(data,mapping=aes(x=avg_log2FC, y=-log10(p_val_adj),
color=DEG))+
geom_vline(xintercept=c(-1,0,1),color="grey",linetype=2)+
geom_hline(yintercept=-log10(0.05),color="grey",linetype=2)+
xlab("log2(FC)")+ylab("-log10(FDR)")+
theme_bw()+
guides(color=guide_legend(title=var1))+
#ggtitle("Low VS High HAPS")+
 theme(legend.position=c(0.5,0.5),
text=element_text(size=15),
plot.title=element_text(size=15,hjust=0.5),
#panel.border=element_blank(),
axis.line.x=element_line(color="black"))+
geom_point(aes(x=data_posi$avg_log2FC,
y=-log10(data_posi$p_val_adj)),color="pink")+
geom_point(aes(x=data_nege$avg_log2FC,
y=-log10(data_nege$p_val_adj)),color="lightblue")+
geom_point(aes(x=data_un$avg_log2FC,
y=-log10(data_un$p_val_adj)),color="grey")+
scale_color_manual(values=c("lightblue","grey","pink"),
labels=c(paste0("down: ",down_num),
paste0("unchange: ",unchange_num),paste0("up: ",up_num)))+

geom_label_repel(data=data_posi,aes(x=avg_log2FC, 
y=-log10(data_posi$p_val_adj),
label =gene),
point.padding=unit(1, "lines"),
max.overlaps=100,
segment.colour = "grey",fill="pink", size=3,
fontface="bold",,color="black",
arrow = arrow(length=unit(0.01,"npc")))+

geom_label_repel(data=data_nege,aes(x=avg_log2FC, 
y=-log10(data_nege$p_val_adj),
label =gene),
point.padding=unit(1, "lines"),
max.overlaps=100,
segment.colour = "grey",fill="lightblue",size=3,
fontface="bold",color="black",
arrow = arrow(length=unit(0.01,"npc")))


}


plt<-plot_grid(li[[1]],li[[2]],li[[3]],li[[4]],ncol=1)

pdf("Volcano_0.05_1.pdf",height=10,width=4)
plt
dev.off()

