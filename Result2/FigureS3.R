setwd("D:/Microbiome/Disease_micorbiome_pairs/Merged")

all<-read.xlsx("../DataSheet.xlsx",sheet="all",startRow = 2,sep.names=" ")
dim(all)
# 993 136

all[1:4,1:4]


data<-read.xlsx("../DataSheet.xlsx",sheet="pheno")
data1<-data[,c("ID","age")]

all<-merge(data1,all,by="ID")
dim(all)

colnames(all)[3:45]
colnames(all)[46:137]

need<-all[,c("age",colnames(all)[46:137])]
need[1:4,1:4]
need<-na.omit(need)
dim(need)
# 957 93
need[1:4,1:4]
write.xlsx(need,"Data_Age_Microbiome.xlsx",overwrite=TRUE)
bak<-need
#####################


need<-bak
dim(need)

need$age<-paste0("Age",need$age)
table(need$age)

tapply(need$Human_endogenous_retrovirus_K,need$age,summary)

compare_means(Human_endogenous_retrovirus_K~age,data=need)

compaired<-list(c("Age1","Age2"),c("Age1","Age3"),
c("Age1","Age4"),c("Age1","Age5"),
c("Age2","Age3"),c("Age2","Age4"),c("Age2","Age5"),
c("Age3","Age4"),c("Age3","Age5"),c("Age4","Age5"))


p1<-ggplot(need,aes(x=age,y=Human_endogenous_retrovirus_K,
fill=age))+
geom_boxplot(width=0.5)+
labs(x='', y= '',title='Human endogenous retrovirus K')+
#geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
scale_fill_tableau()+
theme_bw()+
theme(legend.position='none',
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))+
scale_x_discrete(limits=names(table(need$age)),
labels=c(paste0("Age1 (20-29)","\n",table(need$age)[[1]]),
paste0("Age2 (30-34.9)","\n",table(need$age)[[2]]),
paste0("Age3 (35-39.9)","\n",table(need$age)[[3]]),
paste0("Age4 (40-44.9)","\n",table(need$age)[[4]]),
paste0("Age5 (45+)","\n",table(need$age)[[5]])
))
ggsave("AgeGroup_HERVK.pdf",p1,height=5,width=7)
ggsave("AgeGroup_HERVK.png",p1,height=5,width=7)

######################

need<-bak
dim(need)
# 957 93
need[1:4,1:4]

H<-data.frame()
li<-list()
for(i in 2:dim(need)[2]){
x<-need$age
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
write.csv(H,"Linear_Regression_Age.csv",row.names=F)
sig<-H[H$p<0.05,]
dim(sig)
head(sig)
colnames(sig)<-c("Index","Taxa","Statistic","Pvalue")
write.csv(sig,"Linear_Regression_Age_Sig.csv",row.names=F)
sig<-read.csv("Linear_Regression_Age_Sig.csv")
sig$logp<- -log10(sig$Pvalue)
head(sig)
sig<-sig[order(sig$Statistic),]
label<-sig$Taxa

p<-ggplot(sig,aes(x=Statistic,y=Taxa,fill=logp))+
geom_bar(stat="identity",position="dodge")+
theme_bw()+
scale_fill_gradient(low = "blue",high = "red",na.value = "white")+
scale_y_discrete(limits=label,labels=label)+
ylab("")+xlab("Statistic")+
ggtitle("Linear Regression between Age and Species")+
theme(legend.position="none",
plot.title = element_text(hjust = 0.5),
#axis.text.y=element_text(hjust=1,size=5),
text=element_text(hjust = 0.5))
ggsave("Linear_Regression_Age_plt.pdf",p,height=7,width=7)
ggsave("Linear_Regression_Age_plt.png",p,height=7,width=7)


for(j in 1:dim(sig)[1]){
i<-sig[j,1]
x<-need$age
y <- need[,i]
relation <- lm(y~x)
summary(summary(relation))
t<-summary(relation)$coefficients[2,3]
p<-summary(relation)$coefficients[2,4]
#png(paste0(colnames(need)[i],".png"),height=500,width=500)
pdf(paste0(colnames(need)[i],".pdf"),height=5,width=5)
plot(x,y,col = "blue",
main = paste0("Linear Regression","\n",
"t value = ",round(t,3),"\n",
"Pr(>|t|) = ",round(p,3)),
abline(lm(y~x)),cex = 1.3,pch = 16,
xlab = "Age",ylab = colnames(need)[i])
dev.off()
}


###################################

need<-bak
need$AgeGroup<-rep("",dim(need)[1])
need[which(need$age==1),]$AgeGroup<-"Age [20-30)"
need[which(need$age!=1),]$AgeGroup<-"Age [30+)"
table(need$AgeGroup)
need1<-melt(need,id.vars=c("age","AgeGroup"))
length(table(need1$variable))

stat<-compare_means(value~AgeGroup,need1,group.by="variable")
table(stat$p.signif)
write.csv(stat,"Wilcox_AgeGroup_30.csv",row.names=F)
taxa<-data.frame(variable=stat[stat$p.signif!="ns",]$variable)
need2<-merge(taxa,need1,by="variable")

my_comparisons<-list(c("Age [20-30)","Age [30+)"))
plt<-ggplot(need2,aes(x=AgeGroup,y=value,color=AgeGroup))+
geom_boxplot()+geom_point()+
scale_colour_nejm()+
theme_bw()+
theme(legend.position="none",text=element_text(size=15))+
facet_wrap(.~variable,scales="free_y")+
stat_compare_means(comparisons = my_comparisons)
h<-ceiling((dim(taxa)[1])/4)
ggsave("Wilcox_AgeGroup_30.pdf",plt,height=3*h,width=12)


###############


need<-bak
need$AgeGroup<-rep("",dim(need)[1])
need[which(need$age==1),]$AgeGroup<-"Age [20-30)"
need[which(need$age!=1),]$AgeGroup<-"Age [30+)"
table(need$AgeGroup)
need1<-melt(need,id.vars=c("age","AgeGroup"))
length(table(need1$variable))
stat<-compare_means(value~AgeGroup,need1,group.by="variable")
table(stat$p.signif)
write.csv(stat,"Wilcox_AgeGroup_30.csv",row.names=F)
taxa<-data.frame(variable=stat[stat$p.signif!="ns",]$variable)
need2<-merge(taxa,need1,by="variable")
my_comparisons<-list(c("Age [20-30)","Age [30+)"))
label<-sort(unique(as.vector(need2$variable)),decreasing=T)
plt1<-ggplot(need2,aes(x=variable,y=value,fill=AgeGroup))+
geom_boxplot(width=0.4)+
labs(x='', y= 'Percentage',title='')+
scale_fill_nejm()+
theme_bw()+
stat_compare_means(aes(group = AgeGroup),
method = "wilcox.test",
label = "p.signif",color="red",label.y =42,
symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
symbols = c("***", "**", "*", "ns")))+
ylim(c(0,45))+
geom_signif(comparisons = my_comparisons,
step_increase = 0.1,map_signif_level = T,test = wilcox.test)+
#coord_flip()+
scale_x_discrete(limits=label,labels=label)+
theme(legend.position=c(0.8,0.7),
text=element_text(size=15),
axis.text.x=element_text(vjust=0.5,angle=45),
plot.title = element_text(size=15,hjust = 0.5))
ggsave("Wilcox_AgeGroup_30_variable.pdf",plt1,height=5,width=9)

###################################

need<-bak
need$AgeGroup<-rep("",dim(need)[1])
need[which(need$age==1),]$AgeGroup<-"Age [20-30)"
need[which(need$age==2 | need$age==3),]$AgeGroup<-"Age [30-40)"
need[which(need$age==4 | need$age==5),]$AgeGroup<-"Age [40+)"
table(need$AgeGroup)
need1<-melt(need,id.vars=c("age","AgeGroup"))
length(table(need1$variable))

stat<-compare_means(value~AgeGroup,need1,group.by="variable")
table(stat$p.signif)
write.csv(stat,"Wilcox_AgeGroup_30_40.csv",row.names=F)
taxa<-data.frame(variable=unique(stat[stat$p.signif!="ns",]$variable))
need2<-merge(taxa,need1,by="variable")

my_comparisons<-list(c("Age [20-30)","Age [30-40)"),
c("Age [30-40)","Age [40+)"),
c("Age [20-30)","Age [40+)"))

plt<-ggplot(need2,aes(x=AgeGroup,y=value,color=AgeGroup))+
geom_boxplot()+geom_point()+
scale_colour_nejm()+
theme_bw()+
theme(legend.position="none",text=element_text(size=15))+
facet_wrap(.~variable,scales="free_y")+
stat_compare_means(comparisons = my_comparisons)
h<-ceiling((dim(taxa)[1])/4)
ggsave("Wilcox_AgeGroup_30_40.pdf",plt,height=3*h,width=14)

####################################

need<-bak
need$AgeGroup<-rep("",dim(need)[1])
need[which(need$age==1),]$AgeGroup<-"Age [20-30)"
need[which(need$age==2 | need$age==3),]$AgeGroup<-"Age [30-40)"
need[which(need$age==4 | need$age==5),]$AgeGroup<-"Age [40+)"
table(need$AgeGroup)
need1<-melt(need,id.vars=c("age","AgeGroup"))
length(table(need1$variable))
stat<-compare_means(value~AgeGroup,need1,group.by="variable")
table(stat$p.signif)
write.csv(stat,"Wilcox_AgeGroup_30_40.csv",row.names=F)
taxa<-data.frame(variable=unique(stat[stat$p.signif!="ns",]$variable))
need2<-merge(taxa,need1,by="variable")
stat<-compare_means(value~AgeGroup,need2,group.by="variable")
table(stat$p.signif)
write.csv(stat,"stat.csv")
stat<-read.csv("stat.csv")
stat$com<-paste0(stat$group1,"_",stat$group2)
format<-dcast(stat,variable~com,value.var="p.signif")
label<-as.vector(format$variable)
lab<-data.frame(Inx=1:length(label),variable=label)
need<-merge(lab,format,by="variable")
need<-need[order(need$Inx),c(1,2,4,3,5)]
need[is.na(need)]<-"ns"
text.format<-paste0(need[,3],"\n",need[,4],"\n",need[,5])
my_comparisons<-list(c("Age [20-30)","Age [30-40)"),
c("Age [20-30)","Age [40+)"),
c("Age [30-40)","Age [40+)"))
plt2<-ggplot(need2,aes(x=variable,y=value,color=AgeGroup))+
geom_boxplot(outlier.size=0.01)+
xlab("")+ylab("")+theme_bw()+
ylim(c(0,43))+
stat_compare_means(comparisons = my_comparisons)+
scale_x_discrete(limits=label,labels=label)+
theme(legend.position=c(0.6,0.7),
text=element_text(hjust = 0.5),
axis.text.x=element_text(angle=45,hjust=1),
plot.title=element_text(hjust = 0.5))+
# scale_fill_manual(values = cols2)+
scale_colour_nejm()+
annotate("text", x = c(1:12), y = 40, size = 4, 
colour = "red",label=text.format)
ggsave("Wilcox_AgeGroup_30_40_variable.pdf",plt2,height=5,width=9)



###################################


need<-bak
need$AgeGroup<-rep("",dim(need)[1])
need[which(need$age==1 | need$age==2 | need$age==3),]$AgeGroup<-"Age [20-40)"
need[which(need$age==4 | need$age==5),]$AgeGroup<-"Age [40+)"
table(need$AgeGroup)
need1<-melt(need,id.vars=c("age","AgeGroup"))
length(table(need1$variable))

stat<-compare_means(value~AgeGroup,need1,group.by="variable")
table(stat$p.signif)
write.csv(stat,"Wilcox_AgeGroup_40.csv",row.names=F)
taxa<-data.frame(variable=unique(stat[stat$p.signif!="ns",]$variable))
need2<-merge(taxa,need1,by="variable")

my_comparisons<-list(c("Age [20-40)","Age [40+)"))

plt<-ggplot(need2,aes(x=AgeGroup,y=value,color=AgeGroup))+
geom_boxplot()+geom_point()+
scale_colour_nejm()+
theme_bw()+
theme(legend.position="none",text=element_text(size=15))+
facet_wrap(.~variable,scales="free_y")+
stat_compare_means(comparisons = my_comparisons)
h<-ceiling((dim(taxa)[1])/4)
ggsave("Wilcox_AgeGroup_40.pdf",plt,height=3*h,width=4)

###################################
