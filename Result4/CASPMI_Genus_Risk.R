setwd("D:/Microbiome/5.TCGA")

match<-read.xlsx("Shogun_Filter_mpa_table.xlsx",sheet=2)
dim(match)
head(match)
match<-match[,c("Taxa","Genus")]

data<-read.xlsx("Shogun_Filter_mpa_table.xlsx",sheet=1)
dim(data)
data[1:6,1:6]

all<-merge(match,data,by="Taxa")
dim(all);dim(data)
all[1:6,1:6]

var<-as.data.frame(table(all$Genus))
dim(var)
head(var)

i=1
var1<-as.vector(var[i,1])
all1<-all[all$Genus==var1,-c(1,2)]
all2<-t(data.frame(apply(all1,2,mean)))
all3<-data.frame(Genus=var1,all2)
H<-all3
for(i in 2:dim(var)[1]){
var1<-as.vector(var[i,1])
all1<-all[all$Genus==var1,-c(1,2)]
all2<-t(data.frame(apply(all1,2,mean)))
all3<-data.frame(Genus=var1,all2)
H<-rbind(H,all3)
}
write.csv(H,"Genus_samples.csv",row.names=F)

H<-read.csv("Genus_samples.csv",row.names=1)
dim(H)
H[1:4,1:4]
colnames(H)<-gsub("X","",colnames(H))
mpa<-data.frame(t(H))
mpa[1:4,1:4]
write.csv(taxa,"Genus_taxa.csv")
mpa$Sample<-rownames(mpa)

group<-read.table("Group",header=T,sep="\t")
head(group)
al<-merge(group,mpa,by="Sample")
al[1:10,1:10]
write.csv(al,"Genus_Group.csv",row.names=F)
data1<-melt(al,id.vars=1:4)
head(data1)

var<-as.data.frame(table(data1$variable))
H<-data.frame()
for(i in 1:dim(var)[1]){
var1<-as.vector(var[i,1])
data2<-data1[as.vector(data1$variable==var1),]
g<-names(sort(tapply(data2$value,data2$Group,sum),decreasing=T))[1]
g1<-cbind(var1,g)
H<-rbind(H,g1)
}
write.csv(H,"Genus_Group_sum.csv",row.names=F)

###########

H<-read.csv("Genus_Group_sum.csv")
taxa<-read.csv("Genus_samples.csv",row.names=1)
colnames(taxa)<-gsub("X","",colnames(taxa))
dim(taxa)
taxa[1:4,1:4]

G2_species<-as.vector(H[H[,2]=="G2",1])
G3_species<-as.vector(H[H[,2]=="G3",1])

G2_species<-intersect(G2_species,rownames(taxa))
G3_species<-intersect(G3_species,rownames(taxa))
length(G2_species)
length(G3_species)


G2_metagenome <- taxa[G2_species, ]
G3_metagenome <- taxa[G3_species, ]
dim(G2_metagenome)
dim(G3_metagenome)

G2_metagenome[1:4,1:4]
G3_metagenome[1:4,1:4]

G3_shannon<-vegan::diversity(t(G3_metagenome),'shannon')
G2_shannon<-vegan::diversity(t(G2_metagenome),'shannon')


R_G3 <- apply(G3_metagenome, 2, function(i) (sum(i > 0))) 
R_G2 <- apply(G2_metagenome, 2, function(i) (sum(i > 0)))


G3_prime <- length(G3_species)
G2_prime <- length(G2_species)

psi_G3 <- ((R_G3/G3_prime)*G3_shannon) 
psi_G2 <- ((R_G2/G2_prime)*G2_shannon)

risk_score <- data.frame(log10((psi_G2+0.00001)/(psi_G3+0.00001))) # 0.00001 added to avoid having the denominator as 0
colnames(risk_score)<-"score"
risk_score$Sample<-rownames(risk_score)

mer<-data.frame(G3_shannon=G3_shannon,G2_shannon=G2_shannon,
R_G3=R_G3,R_G2=R_G2,
psi_G3=psi_G3,psi_G2=psi_G2,risk_score=risk_score)
write.csv(mer,"CASMPI_Risk_Score.csv",row.names=F)

res<-merge(group,risk_score,by="Sample")
dim(group);dim(risk_score);dim(res)

res$new<- rep("Low Risk",dim(res)[1])
res[(res$score>0),]$new<-"High Risk"
table(res$new)
st<-table(res$new,res$Group)
p<-fisher.test(st)$p.value
st1<-as.data.frame(st)
colnames(st1)<-c("Risk","Group","Freq")
write.csv(res, "Risk_score_Group.csv",row.names=F)


score<-res$score
pdf("Density_Risk_Group.pdf",height=4,width=4)
r <- hist(score, col = "lightblue", border = "grey50",
freq = F,ylim=c(0,1))
text(r$mids, r$density, r$counts, pos=3,col = "black")
lines(density(score), col = "red")
dev.off()

summary(res$score)


count_table<-c("[-5,-4)","[-4,-3)","[-3,-2)",
"[-2,-1)","[-1,0)","[0,1)",
"[1,2)","[2,3)","[3,4)","[4,5)")

res$bins<-rep("",dim(res)[1])

res[which(res$score>= -5 & res$score< -4),]$bins<-count_table[1]
res[which(res$score>= -4 & res$score< -3),]$bins<-count_table[2]
res[which(res$score>= -3 & res$score< -2),]$bins<-count_table[3]
res[which(res$score>= -2 & res$score< -1),]$bins<-count_table[4]
res[which(res$score>= -1 & res$score< 0),]$bins<-count_table[5]
res[which(res$score>= 0 & res$score< 1),]$bins<-count_table[6]
res[which(res$score>= 1 & res$score< 2),]$bins<-count_table[7]
res[which(res$score>= 2 & res$score< 3),]$bins<-count_table[8]
res[which(res$score>= 3 & res$score< 4),]$bins<-count_table[9]
res[which(res$score>= 4 & res$score< 5),]$bins<-count_table[10]

st<-table(res$bins,res$Group)
write.csv(st,"res_risk_score_bin.csv")
st<-read.csv("res_risk_score_bin.csv",row.names=1)

st[1,]<-st[1,]/sum(st[1,])
st[2,]<-st[2,]/sum(st[2,])
st[3,]<-st[3,]/sum(st[3,])
st[4,]<-st[4,]/sum(st[4,])
st[5,]<-st[5,]/sum(st[5,])
st[6,]<-st[6,]/sum(st[6,])
st[7,]<-st[7,]/sum(st[7,])
st[8,]<-st[8,]/sum(st[8,])
write.csv(st,"res_risk_score_bin1.csv")


#######

pdf("Risk_bins.pdf",height=6,width=6)
hist(score,xaxt='n',yaxt='n',main = "",
axes=F, xlab="Risk Score", ylab=NA)
axis(side = 4)
axis(side = 1)

par(new=T)
plot(st$G1, type="o", col="#F8766D", pch="o", lty=1,
xaxt='n',yaxt='n',xlab = NA,ylab = "Proportion of Samples (%)",
ylim=c(0,1))
par(new=T)
plot(st$G2, type="o", col="#7CAE00", pch="o", lty=1,
xaxt='n',yaxt='n',xlab = NA,ylab = "Proportion of Samples (%)",
ylim=c(0,1))
par(new=T)
plot(st$G3, type="o", col="#00BFC4", pch="o", lty=1,
xaxt='n',yaxt='n',xlab =NA,ylab = "Proportion of Samples (%)",
ylim=c(0,1))
axis(side=2)
par(new=T)
plot(st$G4, type="o", col="#C77CFF", pch="o", lty=1,
xaxt='n',yaxt='n',xlab =NA,ylab = "Proportion of Samples (%)",
ylim=c(0,1))
axis(side=2)
dev.off()

##############


plt1<-ggplot(res,aes(x=score))+
geom_density(alpha=0.8)+
theme_bw()+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")

plt2<-ggplot(res, aes(x = pc_x, y = pc_y,col=new))+
geom_point(alpha=0.8)+
theme_bw()+
stat_ellipse(level = 0.95)+
scale_colour_manual(values = c("#CB3425","#3F5688"))+
geom_vline(xintercept=0,linetype=2,color="grey")+
geom_hline(yintercept=0,linetype=2,color="grey")+
theme(
plot.title = element_text(size=15,hjust = 0.5),
legend.position='right')+
labs(color="Risk Group")

plt<-plot_grid(plt1,plt2,rel_widths=c(1,1.3))

ggsave("Risk_Group_density.png",plt,height=5,width=10)
ggsave("Risk_Group_density.pdf",plt,height=5,width=10)


##################


compaired<-list(c("G1","G2"),c("G2","G3"),c("G3","G4"),
c("G1","G3"),c("G2","G4"),
c("G1","G4"))


p1<-ggplot(res,aes(x=Group,y=score,fill=Group))+
geom_boxplot(width=0.4)+
labs(x='Microbiome Group', y= "Risk score")+
geom_signif(comparisons = compaired,
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
theme_bw()+
theme(legend.title=element_blank(),
legend.position='top',
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))


p2<-ggplot(st1,aes(x=Risk,y=Freq,fill=Group))+
geom_bar(stat="identity",position=position_dodge(0.8),width=0.4)+
theme_bw()+xlab("Risk Group")+ylab("Numbers of Samples")+
annotate("text", x = 1, y = 500, colour = "black",
label=paste0("Fisher's exact test","\n","Pvalue = ",signif(p,3)) )+
theme(legend.position='top',
legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))

p3<-ggplot(res,aes(x=score,fill=new))+
geom_density(alpha=0.8)+
theme_bw()+labs(x='Risk Score', y= "Density of Risk Score")+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position='top',
plot.title = element_text(size=15,hjust = 0.5))

p4<-ggplot(res,aes(x=new,y=score,fill=new))+
geom_boxplot(width=0.4)+
labs(x='Risk Group', y= "Risk score")+
geom_signif(comparisons = list(c("High Risk","Low Risk")),
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
theme_bw()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
legend.position='top',
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))


plt_all<-plot_grid(p1,p2,p3,p4,nrow=1)
plt_all

ggsave("CASPMI_Risk_Group.pdf",plt_all,height=4,width=14)
ggsave("CASPMI_Risk_Group.png",plt_all,height=4,width=14)

################


