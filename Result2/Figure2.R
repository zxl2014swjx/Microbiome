setwd("D:/Microbiome/Combine/Analysis")


data1<-read.table("Phenotype_Group.xls",header=T,sep="\t")

library(reshape2)
data2<-data1[,c(1:6,9:23)]
data2.2<-melt(data2,id.vars=c(1:6))
head(data2.2)
data2.2<-data2.2[which(data2.2$value!="NaN"),]
data2.2$value<-as.numeric(as.vector(data2.2$value))

data2.2.2<-data2.2[which(data2.2$Group=="G2" |
data2.2$Group=="G3"),]
stat<-compare_means(value~Group,data2.2.2,group.by="variable")
compaired1<-list(c("G2","G3"))

var<-as.data.frame(stat[stat$p.signif!="ns",]$variable)

li<-list()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data2.3<-data2.2.2[data2.2.2$variable==var1,]
data2.3$label<-rep("NA",dim(data2.3)[1])
num<-table(data2.3$Group)
lim<-names(num)
lab<-paste0(names(num),"\n",as.numeric(as.vector(num)))
li[[i]]<-ggplot(data2.3,aes(x=Group,y=value,fill=Group))+
geom_violin(width=0.8)+
geom_boxplot(width=0.4)+
labs(x='', y= '',title=toupper(var1))+
geom_signif(comparisons = compaired1,
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
scale_fill_tableau()+
scale_x_discrete(limits=lim,labels=lab)+
theme_bw()+
theme(legend.position='none',
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))
}
length(li)

plt<-gridExtra::arrangeGrob(grobs=li,nrow=1)
dev.off()

ggsave("Group_variable_G23.png",plt,height=4,width=14)
ggsave("Group_variable_G23.pdf",plt,height=4,width=14)



##################

data3<-data1[,c(1:6,24:47)]
data3.2<-melt(data3,id.vars=c(1:6))
head(data3.2)
data3.2<-data3.2[which(data3.2$value!="NaN"),]

data3.2.3<-data3.2[which(data3.2$Group=="G2" |
data3.2$Group=="G3"),]

var<-as.data.frame(table(data3.2.3$variable))
H<-data.frame()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data3.3<-data3.2.3[data3.2.3$variable==var1,]
data3.3$variable<-as.vector(data3.3$variable)
st<-table(data3.3$Group,data3.3$value)
p.value<-chisq.test(st)$p.value
c<-t(data.frame(c(var1,p.value)))
H<-rbind(H,c)
}
rownames(H)<-1:dim(H)[1]
colnames(H)<-c('variable','p.value')
var<-H[H$p.value<0.05,]

li<-list()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data3.3<-data3.2.3[data3.2.3$variable==var1,]
data3.3$variable<-as.vector(data3.3$variable)
st<-table(data3.3$Group,data3.3$value)
p.value<-chisq.test(st)$p.value
lim<-rownames(st)
if(dim(st)[2]==2){
lab<-paste0(rownames(st),"\n",st[,1],"\n",st[,2])
}
st[1,]<-st[1,]/sum(st[1,])
st[2,]<-st[2,]/sum(st[2,])
st1<-as.data.frame(st)
colnames(st1)<-c("Group","var1","value")
li[[i]]<-ggplot(st1,aes(x=Group,y=value,fill=var1))+
geom_bar(stat="identity",position="stack")+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
xlab("")+
theme_bw()+labs(y=paste0("Percentage of ",var1),
title=paste0("Fisher p = ",round(p.value,3)))+
#scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(value,3)*100,"%")),
hjust = "center",size=6,color="white")+
#scale_x_discrete(limits=lim,labels=lab)+
#geom_vline(xintercept=c(558),linetype=1,color="black",size=0.01)
theme(text=element_text(size=12,hjust = 0.5),
legend.position="none",
plot.title = element_text(size=15,hjust = 0.5))+
guides(fill = guide_legend(title = ""))
}
multi<-li
length(li)


plt2<-plot_grid(li[[1]],li[[2]],li[[3]],nrow=1)
dev.off()
ggsave("Group_variable_chisq_G23.png",plt2,height=4,width=12)
ggsave("Group_variable_chisq_G23.pdf",plt2,height=4,width=12)


##############
setwd("D:/Microbiome/Combine/Analysis/Species")

data_need<-read.csv('data_need.csv')


data1<-melt(data_need,id.vars=c(1:4))
head(data1)
head(data1)
data1<-data1[data1$value!=0,]

data1<-data1[which(data1$Group=="G2" |
data1$Group=="G3"),]

stat<-compare_means(value~Group,data1,group.by="variable")
dim(stat)
head(stat)

write.xlsx(stat,"Statistic_Species_G23.xlsx")
stat1<-as.data.frame(stat)
head(stat1)
stat1<-stat1[stat1$p.signif!="ns",]
dim(stat1)

compaired<-list(c("G2","G3"))

li<-list()
var<-as.data.frame(table(stat1$variable))
var<-var[var[,2]>0,]
dim(var)
data1$Group<-as.factor(data1$Group)
head(data1$Group)


data2<-merge(data1,stat1,by="variable")

compare_means(value~Group,data2,group.by="variable")

label<-sort(unique(as.vector(data2$variable)),decreasing=T)
data2$value<-data2$value*100

p3.1<-ggplot(data2,aes(x=variable,y=value,fill=Group))+
geom_boxplot(width=0.4)+
labs(x='', y= 'Percentage',title='')+
#scale_fill_manual(values = c("#CB3425","#3F5688"))+
scale_fill_tableau()+
theme_bw()+
stat_compare_means(aes(group = Group),
method = "wilcox.test",
label = "p.signif",color="red",label.y =40,
symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
symbols = c("***", "**", "*", "ns")))+
ylim(c(0,45))+
geom_signif(comparisons = compaired1,
step_increase = 0.1,map_signif_level = T,test = wilcox.test)+
#coord_flip()+scale_x_discrete(limits=label,labels=label)+
theme(legend.position='right',
axis.text.x=element_text(vjust=0.5,angle=90),
plot.title = element_text(size=15,hjust = 0.5))

ggsave("Species_Wilcox_G23.pdf",p3.1,height=6,width=8)
ggsave("Species_Wilcox_G23.png",p3.1,height=6,width=8)







