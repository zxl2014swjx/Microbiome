#Combine Kraken and Shogun Pipeline
cd /p300s/jiapl_group/zhuxl/Microbiome/Jieguo/Result
ls >../Sample1
cd /p300s/jiapl_group/zhuxl/Microbiome/Analysis
ln -s ../Jieguo/Sample1 ./

R
################

setwd("/p300s/jiapl_group/zhuxl/Microbiome/Analysis")
data<-read.table("Sample1",header=F,sep="\t")
#[1] 988   1
dim(data)
head(data)


var1<-c("kingdom","phylum","class","order","family","genus","species")
for(j in 1:length(var1)){
H1<-data.frame()
var<-as.vector(data[1,1])
a<-read.table(paste0("Result/",var,"/redist.",as.vector(var1[j]),".txt.count"),header=F,sep="\t")
a[,2]<-a[,2]/sum(a[,2])
colnames(a)<-c("Taxa",var)
stat<-data.frame(Sample=var,stat=dim(a)[1])
H1<-rbind(H1,stat)
H<-a
for(i in 2:dim(data)[1]){
var<-as.vector(data[i,1])
a<-read.table(paste0("Result/",var,"/redist.",as.vector(var1[j]),".txt.count"),header=F,sep="\t")
a[,2]<-a[,2]/sum(a[,2])
colnames(a)<-c("Taxa",var)
stat<-data.frame(Sample=var,stat=dim(a)[1])
H<-merge(H,a,all=T)
H1<-rbind(H1,stat)
#print(i)
}
H[is.na(H)]<-0
print(paste0(as.vector(var1[j])," SUCCESS"))
write.table(H,paste0("Shogun/shogun.",as.vector(var1[j]),".table"),row.names=F,sep="\t",quote=F)
write.table(H1,paste0("Shogun/shogun.",as.vector(var1[j]),".stat"),row.names=F,sep="\t",quote=F)
}

########################
data1<-read.table("Shogun/shogun.genus.table",header=T,sep="\t")
data1<-data1[-grep("BacteriaPlasmid",data1$Taxa),]
colnames(data1)<-gsub("^X","",colnames(data1))
dim(data1)
head(data1[,1:4])
st1<-c()
for(i in 1:dim(data1)[1]){
st<-unlist(strsplit(as.vector(data1[i,1]),";"))[6]
st1<-c(st1,st)
}
st1<-gsub("g__","",st1)
rownames(data1)<-as.vector(data1$Taxa)
data1$Taxa<-st1
data1[1:4,1:4]


data2<-read.table("Kraken/kraken.genus.table",header=T,sep="\t")
colnames(data2)<-gsub("^X","",colnames(data2))
data2[1:4,1:4]
data2$Taxa<-as.vector(data2$Taxa)
data2$Taxa<-gsub(" ","_",data2$Taxa)
dim(data2)
head(data2[,1:4])
for(j in 2:dim(data2)[2]){
data2[,j]<-data2[,j]/sum(data2[,j])
}
data2[1:4,1:4]

need<-intersect(colnames(data1),colnames(data2))
data1<-data1[,need]
data2<-data2[,need]
table(colnames(data1)==colnames(data2))



a<-data1[,c(1,2)]
b<-data2[,c(1,2)]
c<-merge(a,b,by="Taxa")
#d<-data.frame(Taxa=c[,1],stat=(c[,2]+c[,3])/2)
d<-data.frame(Taxa=c[,1],stat=c[,3])
colnames(d)[2]<-as.vector(colnames(data1)[2])
H<-d

for(k in 3:dim(data1)[2]){
a<-data1[,c(1,k)]
b<-data2[,c(1,k)]
c<-merge(a,b,by="Taxa")
d<-data.frame(Taxa=c[,1],stat=c[,3])
colnames(d)[2]<-as.vector(colnames(data1)[k])
if(d[,1]==H[,1]){
d1<-data.frame(d[,-1])
colnames(d1)<-colnames(d)[2]
H<-cbind(H,d1)
}}
H[1:4,1:4]
write.csv(H,"H.csv")


need<-data.frame(rownames(data1),data1$Taxa)
colnames(need)<-c("Taxa1","Taxa")
write.table(need,"Match_Genus.xls",row.names=F,sep="\t",quote=F)


al<-merge(need,H)
al[1:4,1:4]
write.table(al,"Kraken_Shogun_Genus.tsv",sep="\t",quote=F)
al2<-al[,-c(1:2)]
table(rowSums(al2)==0)
rownames(al2)<-as.vector(al$Taxa)
al2[1:4,1:4]
write.table(al2,"Kraken_Genus.tsv",sep="\t",quote=F)
al3<-al2
rownames(al3)<-as.vector(al$Taxa1)
al3[1:4,1:4]
write.table(al3,"Shogun_Genus.tsv",sep="\t",quote=F)


mpa_table<-al2
mpa_table[1:4,1:4]
mpa_table[mpa_table<=0.01]<-0
mpa_table<-mpa_table[rowSums(mpa_table)!=0,colSums(mpa_table)!=0]
dim(mpa_table)
write.table(mpa_table,"Filter_Kraken_Genus.tsv",sep="\t",quote=F)

mpa_table<-al3
mpa_table[1:4,1:4]
mpa_table[mpa_table<=0.01]<-0
mpa_table<-mpa_table[rowSums(mpa_table)!=0,colSums(mpa_table)!=0]
dim(mpa_table)
write.table(mpa_table,"Filter_Shogun_Genus.tsv",sep="\t",quote=F)


########################
data1<-read.table("Shogun/shogun.species.table",header=T,sep="\t")
data1<-data1[-grep("BacteriaPlasmid",data1$Taxa),]
colnames(data1)<-gsub("^X","",colnames(data1))
dim(data1)
head(data1[,1:4])
st1<-c()
for(i in 1:dim(data1)[1]){
st<-unlist(strsplit(as.vector(data1[i,1]),";"))[7]
st1<-c(st1,st)
}
st1<-gsub("s__","",st1)
rownames(data1)<-as.vector(data1$Taxa)
data1$Taxa<-st1
data1[1:4,1:4]



data2<-read.table("Kraken/kraken.species.table",header=T,sep="\t")
colnames(data2)<-gsub("^X","",colnames(data2))
data2[1:4,1:4]
data2$Taxa<-as.vector(data2$Taxa)
data2$Taxa<-gsub(" ","_",data2$Taxa)
dim(data2)
head(data2[,1:4])
for(j in 2:dim(data2)[2]){
data2[,j]<-data2[,j]/sum(data2[,j])
}
data2[1:4,1:4]

need<-intersect(colnames(data1),colnames(data2))
data1<-data1[,need]
data2<-data2[,need]
table(colnames(data1)==colnames(data2))



a<-data1[,c(1,2)]
b<-data2[,c(1,2)]
c<-merge(a,b,by="Taxa")
#d<-data.frame(Taxa=c[,1],stat=(c[,2]+c[,3])/2)
d<-data.frame(Taxa=c[,1],stat=c[,3])
colnames(d)[2]<-as.vector(colnames(data1)[2])
H<-d

for(k in 3:dim(data1)[2]){
a<-data1[,c(1,k)]
b<-data2[,c(1,k)]
c<-merge(a,b,by="Taxa")
d<-data.frame(Taxa=c[,1],stat=c[,3])
colnames(d)[2]<-as.vector(colnames(data1)[k])
if(d[,1]==H[,1]){
d1<-data.frame(d[,-1])
colnames(d1)<-colnames(d)[2]
H<-cbind(H,d1)
}}
H[1:4,1:4]
write.csv(H,"H.csv")

need<-data.frame(rownames(data1),data1$Taxa)
colnames(need)<-c("Taxa1","Taxa")
write.table(need,"Match_Genus.xls",row.names=F,sep="\t",quote=F)



al<-merge(need,H)
al[1:4,1:4]
write.table(al,"Kraken_Shogun_Species.tsv",sep="\t",quote=F)
al2<-al[,-c(1:2)]
table(rowSums(al2)==0)
rownames(al2)<-as.vector(al$Taxa)
al2[1:4,1:4]
write.table(al2,"Kraken_Species.tsv",sep="\t",quote=F)
al3<-al2
rownames(al3)<-as.vector(al$Taxa1)
al3[1:4,1:4]
write.table(al3,"Shogun_Species.tsv",sep="\t",quote=F)


mpa_table<-al2
mpa_table[1:4,1:4]
mpa_table[mpa_table<=0.01]<-0
mpa_table<-mpa_table[rowSums(mpa_table)!=0,colSums(mpa_table)!=0]
dim(mpa_table)
write.table(mpa_table,"Filter_Kraken_Species.tsv",sep="\t",quote=F)

mpa_table<-al3
mpa_table[1:4,1:4]
mpa_table[mpa_table<=0.01]<-0
mpa_table<-mpa_table[rowSums(mpa_table)!=0,colSums(mpa_table)!=0]
dim(mpa_table)
write.table(mpa_table,"Filter_Shogun_Species.tsv",sep="\t",quote=F)

