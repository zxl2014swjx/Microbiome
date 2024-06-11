setwd("D:/Microbiome/Combine")

dir()


data<-read.table("kraken.species.table",header=T,sep="\t",row.names=1)
filter<-read.table("mpa_table.xls",header=T,sep="\t",row.names=1)
colnames(data)<-gsub("X","",colnames(data))
colnames(filter)<-gsub("X","",colnames(filter))

dim(data)
dim(filter)

data[1:4,1:4]
filter[1:4,1:4]

table(colnames(data)==colnames(filter))

H<-data.frame()
for(i in 1:dim(data)[2]){
sam<-colnames(data)[i]
a<-data[,i]
b<-filter[,i]
a1<-a[a!=0]
b1<-b[b!=0]
st<-cbind(sam,length(a1),length(b1))
H<-rbind(H,st)
}
colnames(H)<-c("SampleID","Kraken2","Filtered")
write.csv(H,"Kraken_Filtered_Stat.csv",row.names=F)

H<-read.csv("Kraken_Filtered_Stat.csv")
summary(H$Kraken2)

H$Suspect_Pollution<-H$Kraken2-H$Filtered
H<-H[order(H$Kraken2,decreasing=T),]
label<-H$SampleID

H1<-melt(H,id.vars=1:2)

plt<-ggplot(H1,aes(x=SampleID,y=value,fill=variable))+
geom_bar(stat="identity",position="stack")+
theme_bw()+labs(y="Percentage")+
scale_y_continuous(expand=c(0,0))+
scale_fill_npg()+
#scale_x_discrete(limits=label,labels=label)+
theme(text=element_text(size=15,hjust = 0.5),
axis.text.x=element_blank(),
axis.ticks.x=element_blank(),
legend.position=c(0.4,0.8))
ggsave("Kraken_Filtered_Stat.png",plt,height=5,width=8)
ggsave("Kraken_Filtered_Stat.pdf",plt,height=5,width=8)





