setwd("D:/Microbiome/3.LEFSE")
library(ggpubr)
library(cowplot)

data<-read.table("Interact_Shogun_mpa_table.res1",header=F,sep="\t")
dim(data)
head(data)

colnames(data)<-c("Biomarker","Abundance","Group","LDA","Pval")
table(data$Group)
data<-data[grep("G",data$Group),]
write.xlsx(data,"LDA_Group.xlsx")
plot(density(data$LDA))

data<-read.xlsx("LDA_Group.xlsx")
table(data$Level,data$Group)



data1<-data[data$Level=="S",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p1<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
#geom_text(aes(label = Species),col="black",vjust = "left")+
theme_bw()+labs(x="LDA Score",y="Species")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="top",
text=element_text(hjust = 0.5))


data1<-data[data$Level=="G",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p2<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
theme_bw()+labs(x="LDA Score",y="Genus")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="none",
text=element_text(hjust = 0.5))

data1<-data[data$Level=="F",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p3<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
theme_bw()+labs(x="LDA Score",y="Family")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="none",
text=element_text(hjust = 0.5))

data1<-data[data$Level=="O",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p4<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
theme_bw()+labs(x="LDA Score",y="Order")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="none",
text=element_text(hjust = 0.5))


data1<-data[data$Level=="C",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p5<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
theme_bw()+labs(x="LDA Score",y="Class")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="none",
text=element_text(hjust = 0.5))

data1<-data[data$Level=="P",]
need<- data1 %>% group_by(Group) %>% top_n(n = 5, wt = LDA)
table(need$Group)
need<-need[order(need$LDA),]
need<-need[order(need$Group,decreasing=T),]
label<-need$label
p6<-ggplot(need,aes(x=LDA,y=label,fill=Group))+
geom_bar(stat="identity",position="dodge",width=0.4)+
scale_y_discrete(limits=label,labels=label)+
theme_bw()+labs(x="LDA Score",y="Phylum")+
scale_fill_manual(values=c("G1"="#F8766D",
"G2"="#7CAE00","G3"="#00BFC4","G4"="#C77CFF"))+
scale_x_continuous(expand=c(0,0))+
theme(legend.position="none",
text=element_text(hjust = 0.5))

plt_all<-plot_grid(p1,p2,p3,p4,p5,p6,
nrow=6,ncol=1,rel_heights=c(2,1.5,1.5,1,1,1))
plt_all

ggsave("Dominant_Taxanomy.pdf",plt_all,height=11.6,width=8.2)
ggsave("Dominant_Taxanomy.png",plt_all,height=11.6,width=8.2)



