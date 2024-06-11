setwd("D:/Microbiome/Disease_micorbiome_pairs/Merged")

all<-read.xlsx("../DataSheet.xlsx",sheet="all",startRow = 2,sep.names=" ")
dim(all)
# 993 136

all[1:4,1:4]

colnames(all)[2:44]
colnames(all)[45:136]


var<-colnames(all)[2:44]
H1<-data.frame()
for(k in 1:length(var)){
var1<-var[k]
need<-all[,c(var1,colnames(all)[45:136])]
need<-na.omit(need)
dim(need)
need[1:4,1:4]
H<-data.frame()
for(i in 2:dim(need)[2]){
x<-need[,1]
y <- need[,i]
relation <- lm(y~x)
summary(summary(relation))
t<-summary(relation)$coefficients[2,3]
p<-summary(relation)$coefficients[2,4]
ne<-cbind(i,colnames(need)[i],t,p)
H<-rbind(H,ne)
}
H<-H[which(H$t!="NaN"),]
H<-H[which(H$p!="NaN"),]
head(H)
dim(H)
write.csv(H,paste0("Linear_Regression_",var1,".csv"),row.names=F)
H<-cbind(var1,H)
H1<-rbind(H1,H)
}
write.table(H1,"Linear_Regression_Continuous.xls",
row.names=F,sep="\t",quote=F)

H1<-read.table("Linear_Regression_Continuous.xls",header=T,sep="\t")
sig<-H1[H1$p<0.05,]
table(sig[,1])
colnames(sig)<-c("Phenotype","Index","Taxa","Statistic","Pvalue")
sig$logp<- -log10(sig$Pvalue)
head(sig)
write.table(sig,"Linear_Regression_Continuous_sig.xls",
row.names=F,sep="\t",quote=F)

sig<-read.table("Linear_Regression_Continuous_sig.xls",header=T,sep="\t")

var<-as.data.frame(table(sig[,1]))
li<-list()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
sig1<-sig[sig$Phenotype==var1,]
sig1<-sig1[order(sig1$Statistic),]
label<-sig1$Taxa
li[[i]]<-ggplot(sig1,aes(x=Statistic,y=Taxa,fill=logp))+
geom_bar(stat="identity",position="dodge")+
theme_bw()+
scale_fill_gradient(low = "blue",high = "red",na.value = "white")+
scale_y_discrete(limits=label,labels=label)+
ylab("")+xlab("Statistic")+
ggtitle(var1)+
theme(legend.position="right",
plot.title = element_text(hjust = 0.5),
text=element_text(hjust = 0.5,size=13))
}

h<-ceiling(length(li)/3)
plt<-gridExtra::arrangeGrob(grobs=li,ncol=3)
ggsave("Linear_Regression_Continuous_sig.pdf",plt,height=3*h,width=14)


###################

H1<-read.table("Linear_Regression_Continuous.xls",header=T,sep="\t")
colnames(H1)<-c("Phenotype","Index","Taxa","Statistic","Pvalue")
H1$format<-rep("",dim(H1)[1])
H1[which(H1$Pvalue>=0.05),]$format<-""
H1[which(H1$Pvalue<0.05 & H1$Pvalue>=0.01),]$format<-"*"
H1[which(H1$Pvalue<0.01 & H1$Pvalue>=0.001),]$format<-"***"
H1[which(H1$Pvalue<0.001),]$format<-"***"
table(H1$format)
write.csv(H1,"Linear_Regression_Continuous_Sig.csv",row.names=F)
H1<-read.csv("Linear_Regression_Continuous_Sig.csv")


pmat<-dcast(H1,Phenotype~Taxa,value.var='Statistic')
rownames(pmat)<-as.vector(pmat[,1])
pmat<-pmat[,-1]
pmat[1:4,1:4]
dim(pmat)

psig<-dcast(H1,Phenotype~Taxa,value.var='format')
rownames(psig)<-as.vector(psig[,1])
psig<-psig[,-1]
psig[1:4,1:4]
dim(psig)


bk1<-seq(min(pmat),0,by=0.01)
bk2<-seq(0,max(pmat),by=0.01)
bk<-c(bk1,bk2)
col_sel<-c(colorRampPalette(colors=c("blue","white"))(which(bk==0)),
colorRampPalette(colors=c("white","red"))(length(bk)-which(bk==0)))


pdf("Linear_Phenotype_Microbiome_pval.pdf",height=14,width=16)
pheatmap(as.matrix(pmat),
display_numbers = as.matrix(psig),
#cluster_rows=F,cluster_cols=F,
cellheight=8,cellwidth=8,fontsize=7,
color=col_sel,breaks=bk)
dev.off()

png("Linear_Phenotype_Microbiome_pval.png",height=1000,width=1200)
pheatmap(as.matrix(pmat),
display_numbers = as.matrix(psig),
#cluster_rows=F,cluster_cols=F,
cellheight=8,cellwidth=8,fontsize=7,
color=col_sel,breaks=bk)
dev.off()

######################

H1<-read.csv("Linear_Regression_Continuous_Sig.csv")
head(H1$Pvalue)
sig<-H1[H1$Pvalue<0.05,]


a1<-as.data.frame(table(sig$Taxa))
colnames(a1)<-c("Taxa","sig_Taxa")
b1<-as.data.frame(table(sig$Phenotype))
colnames(b1)<-c("Phenotype","sig_Phenotype")


Taxa<-data.frame(Taxa=unique(sig$Taxa))
sig1<-merge(Taxa,H1,by="Taxa")

sig1<-merge(sig1,a1,by="Taxa",all.x=T)
sig1<-merge(sig1,b1,by="Phenotype",all.x=T)
sig1[is.na(sig1)]<-0


sig1$col<-paste0(sig1$Taxa," (",sig1$sig_Taxa,")")
sig1$row<-paste0(sig1$Phenotype," (",sig1$sig_Phenotype,")")

table(sig1$format)

write.csv(sig1,"Linear_Regression_Continuous_sig1.csv",row.names=F)
sig1<-read.csv("Linear_Regression_Continuous_sig1.csv")


pmat<-dcast(sig1,row~col,value.var='Statistic')
rownames(pmat)<-as.vector(pmat[,1])
pmat<-pmat[,-1]
pmat[1:4,1:4]
dim(pmat)

psig<-dcast(sig1,row~col,value.var='format')
rownames(psig)<-as.vector(psig[,1])
psig<-psig[,-1]
psig[1:4,1:4]
dim(psig)


bk1<-seq(min(pmat),0,by=0.01)
bk2<-seq(0,max(pmat),by=0.01)
bk<-c(bk1,bk2)
col_sel<-c(colorRampPalette(colors=c("blue","white"))(which(bk==0)),
colorRampPalette(colors=c("white","red"))(length(bk)-which(bk==0)))


pdf("Linear_Phenotype_Microbiome_format.pdf",height=10,width=12)
pheatmap(as.matrix(pmat),
display_numbers = as.matrix(psig),
#cluster_rows=F,cluster_cols=F,
cellheight=8,cellwidth=8,fontsize=7,
color=col_sel,breaks=bk)
dev.off()

png("Linear_Phenotype_Microbiome_format.png",height=600,width=900)
pheatmap(as.matrix(pmat),
display_numbers = as.matrix(psig),
#cluster_rows=F,cluster_cols=F,
cellheight=8,cellwidth=8,fontsize=7,
color=col_sel,breaks=bk)
dev.off()

#############

