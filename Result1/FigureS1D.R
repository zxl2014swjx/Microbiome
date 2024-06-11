setwd("D:/Microbiome/Disease_micorbiome_pairs")

all<-read.table("Phenotype_microbiome_relationship.xls",header=T,sep="\t")
colnames(all)[6:20]

colnames(all)[65:156]

# continous variable: spearman or pearson

H<-data.frame()
for(j in 65:156){
#j=65
a<-table(all[,j]>0)
if(length(a)==2){
pre<-a[2]/993
ma<-data.frame(Taxa=colnames(all)[j],Prevelance=pre)}
if(length(a)==1){
pre<-1
ma<-data.frame(Taxa=colnames(all)[j],Prevelance=pre)}
H<-rbind(H,ma)
}

dim(H)
head(H)

#############
data1<-melt(data,id.vars=1:64)
var<-as.data.frame(table(data1$variable))
class<-data.frame()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data2<-data1[as.vector(data1$variable==var1),]
g<-names(sort(tapply(data2$value,data2$Group,sum),decreasing=T))[1]
g1<-cbind(var1,g)
class<-rbind(class,g1)
}
write.csv(class,"Species_Group_sum.csv",row.names=F)
############

class<-read.csv("Species_Group_sum.csv")
colnames(class)<-c("Taxa","Group")

H<-merge(H,class,by="Taxa")
H<-H[order(H[,2],decreasing=T),]
write.table(H,"Microbiome_Samples_Prevance.xls",row.names=F,sep="\t",quote=F)

label<-as.vector(H[,1])
plt<-ggplot(H,aes(x=Taxa,y=Prevelance,fill=Group))+
geom_bar(stat="identity",position="dodge")+
theme_bw()+
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
legend.position=c(0.8,0.8))+
scale_x_discrete(limits=label)+
xlab("Microbiome")

ggsave("Microbiome_Prevalance.pdf",plt,height=5,width=9)




