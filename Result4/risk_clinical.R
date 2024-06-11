setwd("D:/Microbiome/5.GMHI")

library(openxlsx)
data1<-read.xlsx("Risk_score_Group.xlsx",sheet="Clinical")


st<-table(data1$Group,data1$Risk_Group)
chisq.test(st)
st1<-as.data.frame(st)
colnames(st1)<-c("Group","Risk_Group","Freq")

p.value<-chisq.test(st)$p.value
lim<-rownames(st)
if(dim(st)[2]==2){
lab<-paste0(rownames(st),"\n",st[,1],"\n",st[,2])
}


st[1,]<-st[1,]/sum(st[1,])
st[2,]<-st[2,]/sum(st[2,])
st[3,]<-st[3,]/sum(st[3,])
st[4,]<-st[4,]/sum(st[4,])

st1<-as.data.frame(st)
colnames(st1)<-c("Group","Risk_Group","Freq")

plt1<-ggplot(st1,aes(x=Group,y=Freq,fill=Risk_Group))+
geom_bar(stat="identity",position="fill")+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
xlab("")+
theme_bw()+labs(y=paste0("Percentage of Risk_Group"),
title=paste0("Chisq p = ",format(p.value,digits = 3)))+
#scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Freq,3)*100,"%")),
hjust = "center",size=6,color="white")+
scale_x_discrete(limits=lim,labels=lab)+
theme(text=element_text(size=12,hjust = 0.5),
legend.position="top",
plot.title = element_text(size=15,hjust = 0.5))+
guides(fill = guide_legend(title = ""))

ggsave("Bar_Group Risk_Group.pdf",plt1,height=5,width=5)
ggsave("Bar_Group Risk_Group.png",plt1,height=5,width=5)

######################

data2<-data1[,c(1:4,5:19)]
dim(data2)

li<-list()
for(i in 5:dim(data2)[2]){
need<-data2[,c(2,3,i)]
need<-need[which(need[,3]!="NaN"),]
need[,3]<-as.numeric(as.vector(need[,3]))
colnames(need)<-c("risk_score","risk_Group","value")

need1<-need[need$risk_Group=="High Risk",]
m<-as.numeric(need1$risk_score)
n<-as.numeric(need1$value)
r1<-cor(m,n,method="spearman")
p1<-cor.test(m,n)$p.value

need2<-need[need$risk_Group=="Low Risk",]
m<-as.numeric(need2$risk_score)
n<-as.numeric(need2$value)
r2<-cor(m,n,method="spearman")
p2<-cor.test(m,n)$p.value

loc_x=median(need$risk_score)
loc_y=max(need$value)

m<-as.numeric(need$risk_score)
n<-as.numeric(need$value)
r<-cor(m,n,method="spearman")
p<-cor.test(m,n)$p.value

li[[i]]<-ggplot(need,aes(x=risk_score,y=value,col=risk_Group))+
geom_point(na.rm = TRUE)+
stat_smooth(method =  lm)+
labs(x='Risk Score',y=colnames(data2)[i])+
ggtitle(paste0("ALL r=",round(r,3),"  p=",signif(p,3),"\n",
"High Risk r=",round(r1,3),"  p=",signif(p1,3),"\n",
"Low Risk r=",round(r2,3),"  p=",signif(p2,3)))+
scale_color_manual(values = c("#CB3425","#3F5688"))+
theme_bw()+
theme(legend.position="top",
plot.title = element_text(hjust = 0.5))
}
length(li)
li1<-li[5:19]

plt0<-gridExtra::arrangeGrob(grobs=li1,ncol=5)
dev.off()

ggsave("Risk_Group_spearman_Risk.png",plt0,height=18,width=18)
ggsave("Risk_Group_spearman_Risk.pdf",plt0,height=18,width=18)



#################
library(reshape2)

data2<-data1[,c(1:4,5:19)]
data2.2<-melt(data2,id.vars=c(1:4))
head(data2.2)
data2.2<-data2.2[which(data2.2$value!="NaN"),]
data2.2$value<-as.numeric(as.vector(data2.2$value))

data2.2.2<-data2.2[which(data2.2$Risk_Group=="High Risk" |
data2.2$Risk_Group=="Low Risk"),]
stat<-compare_means(value~Risk_Group,data2.2.2,group.by="variable")
compaired1<-list(c("High Risk","Low Risk"))

var<-as.data.frame(stat[stat$p.signif!="ns",]$variable)
var<-as.data.frame(table(data2.2.2$variable))


li<-list()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data2.3<-data2.2.2[data2.2.2$variable==var1,]
data2.3$label<-rep("NA",dim(data2.3)[1])
num<-table(data2.3$Risk_Group)
lim<-names(num)
lab<-paste0(names(num),"\n",as.numeric(as.vector(num)))
li[[i]]<-ggplot(data2.3,aes(x=Risk_Group,y=value,fill=Risk_Group))+
geom_violin(width=0.8)+
geom_boxplot(width=0.2)+
labs(x='', y= '',title=toupper(var1))+
geom_signif(comparisons = compaired1,map_signif_level = F,test = wilcox.test)+
scale_fill_tableau()+
scale_x_discrete(limits=lim,labels=lab)+
theme_bw()+
theme(legend.position='none',
text=element_text(hjust = 0.5),
plot.title = element_text(hjust = 0.5))
}
length(li)

plt<-gridExtra::arrangeGrob(grobs=li,ncol=5)
dev.off()

ggsave("Risk_Group_variable_High Risk3.png",plt,height=12,width=12)
ggsave("Risk_Group_variable_High Risk3.pdf",plt,height=12,width=12)



##################

data3<-data1[,c(1:4,20:43)]
data3.2<-melt(data3,id.vars=c(1:4))
head(data3.2)
data3.2<-data3.2[which(data3.2$value!="NaN"),]

data3.2.3<-data3.2[which(data3.2$Risk_Group=="High Risk" |
data3.2$Risk_Group=="Low Risk"),]

var<-as.data.frame(table(data3.2.3$variable))


H<-data.frame()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data3.3<-data3.2.3[data3.2.3$variable==var1,]
data3.3$variable<-as.vector(data3.3$variable)
st<-table(data3.3$Risk_Group,data3.3$value)
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
data3.3<-data3.2[data3.2$variable==var1,]
data3.3$variable<-as.vector(data3.3$variable)
st<-table(data3.3$Risk_Group,data3.3$value)
p.value<-chisq.test(st)$p.value
st1<-as.data.frame(st)
colnames(st1)<-c("Group","var1","value")
lim<-rownames(st)
if(dim(st)[2]==2){
lab<-paste0(rownames(st),"\n",st[,1],"\n",st[,2])
}
if(dim(st)[2]>2){
lab<-paste0(rownames(st),"\n",
st[,1],"\n",st[,2],"\n",st[,3],"\n",st[,4],"\n",st[,5])
}
li[[i]]<-ggplot(st1,aes(x=Group,y=value,fill=var1))+
geom_bar(stat="identity",position="fill")+
scale_fill_jco()+xlab("")+
theme_bw()+labs(y="Percentage",
title=paste0(var1,"\n","pvalue = ",round(p.value,3)))+
scale_y_continuous(expand=c(0,0))+
scale_x_discrete(limits=lim,labels=lab)+
#geom_vline(xintercept=c(558),linetype=1,color="black",size=0.01)
theme(text=element_text(size=6,hjust = 0.5),
legend.position="top",
plot.title = element_text(size=15,hjust = 0.5))+
guides(fill = guide_legend(title = ""))
}
multi<-li
length(li)


plt2<-gridExtra::arrangeGrob(grobs=li,ncol=6)
dev.off()
ggsave("Risk_Group_variable_chisq_High Risk3.png",plt2,height=14,width=18)
ggsave("Risk_Group_variable_chisq_High Risk3.pdf",plt2,height=14,width=18)


##############







