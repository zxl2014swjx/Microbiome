setwd("D:/Microbiome/Disease_micorbiome_pairs/Age")

all<-read.table("Phenotype_microbiome_relationship.xls",header=T,sep="\t")
dim(all)

all[1:4,1:4]

table(all$age)
table(all$sex,all$age)

colnames(all)[6:20]
colnames(all)[65:156]

need<-all[,c("age",colnames(all)[65:156])]
need[1:4,1:4]
need<-na.omit(need)
bak<-need

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
