setwd("D:/Microbiome/5.Singapore")

taxa<-read.csv("SG10K_Merged_Taxa_Select.csv",row.names=2)
dim(taxa)
taxa<-taxa[,-1]
taxa[1:4,1:4]
taxa[is.na(taxa)]<-0
taxa<-taxa[,colSums(taxa)!=0]
write.csv(taxa,"taxa_input.csv")


H<-read.csv("../5.GMHI/Species_Group_sum.csv")

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
summary(risk_score)
colnames(risk_score)<-"score"
risk_score$sample<-rownames(risk_score)
head(risk_score)
summary(risk_score$score)


mer<-data.frame(G3_shannon=G3_shannon,G2_shannon=G2_shannon,
R_G3=R_G3,R_G2=R_G2,
psi_G3=psi_G3,psi_G2=psi_G2,risk_score=risk_score)
write.csv(mer,"SG10K_Risk_Score.csv",row.names=F)

###############


setwd("D:/Microbiome/5.Singapore")

taxa<-read.csv("taxa_input.csv",row.names=1)
dim(taxa)
head(taxa)
taxa[1:4,1:4]

mpa_table<-taxa

library(microbiome)
alpha_diversity <- microbiome::alpha(mpa_table, 
index = c("richness_observed","chao1", #Richness
"diversity_shannon","diversity_gini_simpson","dominance_gini","diversity_inverse_simpson","diversity_coverage",
"evenness_camargo","evenness_pielou","evenness_simpson","evenness_evar","evenness_bulla",
"dominance_dbp","dominance_dmn","dominance_absolute","dominance_relative","dominance_simpson","dominance_core_abundance",
"rarity_log_modulo_skewness","rarity_low_abundance","rarity_rare_abundance"))
head(alpha_diversity)
dim(alpha_diversity)

BCI<-t(mpa_table)
BCI[1:4,1:4]

row<-rownames(BCI)
col<-colnames(BCI)
pinpu<-matrix(rep(0,length(row)*length(col)),nrow=length(row))
rownames(pinpu)<-row
colnames(pinpu)<-col
for(i in 1:length(row)){
for(j in 1:length(col)){
pinpu[i,j]=as.numeric(as.vector(BCI[i,j]))
}
}
BCI<-pinpu
write.csv(pinpu,"BCI_pinpu.csv")
BCI<-BCI*100
BCI[1:4,1:4]
write.table(BCI,"Filter_BCI.xls",sep="\t",quote=F)
BCI[1,1]


shannon<-vegan::diversity(BCI,'shannon')
simpson<-vegan::diversity(BCI,'simpson')
invsimpson<-vegan::diversity(BCI,'inv')
unbias.simpson <- vegan::simpson.unb(BCI)
library(vegan)
#alpha <- vegan::fisher.alpha(BCI)
species <- vegan::specnumber(BCI) 
evenness <- species/log(species) 
data<-data.frame(shannon,simpson,invsimpson,
species,evenness)

if(names(table(rownames(data)==rownames(alpha_diversity)))){
al<-cbind(data,alpha_diversity)
}
dim(al)
head(al[,c(1:6)])
ID<-rownames(al)
al<-cbind(ID,al)
write.table(al,"Microbiome_Diversity.xls",sep="\t",quote=F,row.names=F)

###############
data<-read.xlsx("SG10K_Risk_Score.xlsx",sheet=2)
dim(data)
head(data)


shanno<-read.table("Microbiome_Diversity.xls",header=T,sep="\t")
dim(shanno)
head(shanno)
#shanno<-shanno[,c(1,2)]
colnames(shanno)[1]<-c("sample")
all<-merge(shanno,data,by="sample")
dim(all);dim(data)
write.csv(all,"SG10K_shnno_risk.csv",row.names=F)

p1<-ggplot(all,aes(x=risk_score,fill=Group))+
geom_density(alpha=0.8)+
theme_bw()+labs(x='Risk Score', y= "Density of Risk Score")+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position='top',
plot.title = element_text(size=15,hjust = 0.5))


p2<-ggplot(all,aes(x=Group,y=risk_score,fill=Group))+
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

plt<-plot_grid(p1,p2,nrow=1)
plt

ggsave("Risk_plt_SG10K.pdf",plt,height=4,width=8)
ggsave("Risk_plt_SG10K.png",plt,height=4,width=8)


##########

need<-all[,c(1,35,2:27)]
write.csv(need,"Risk_alpha.csv",row.names=F)

need<-read.csv("Risk_alpha.csv")
need1<-melt(need,id.vars=1:2)
head(need1)

stat<-compare_means(value~Group, data=need1,group.by= "variable")
dim(stat)
head(stat)
write.csv(stat,"Risk_alpha_stat.csv",row.names=F)


plt2<-ggplot(need1,aes(x=Group,y=value,fill=Group))+
geom_violin(width=0.7)+
geom_boxplot(width=0.1,fill="white")+
facet_wrap(.~variable,scales="free_y")+
labs(x='Risk Group', y= "Risk score")+
geom_signif(comparisons = list(c("High Risk","Low Risk")),
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
theme_bw()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
legend.position='none',
text=element_text(size=15,hjust = 0.5),
plot.title = element_text(size=15,hjust = 0.5))

ggsave("Risk_plt_SG10K_diversity.pdf",plt2,height=12,width=15)
ggsave("Risk_plt_SG10K_diversity.png",plt2,height=12,width=15)


################

child<-all[all$AgeGroup=="Child",]
m<-as.numeric(child$risk_score)
n<-as.numeric(child$Age)
r1<-cor(m,n,method="spearman")
p1<-cor.test(m,n)$p.value

adult<-all[all$AgeGroup=="Adult",]
m<-as.numeric(adult$risk_score)
n<-as.numeric(adult$Age)
r2<-cor(m,n,method="spearman")
p2<-cor.test(m,n)$p.value

loc_x=median(all$risk_score)
loc_y=max(all$Age)

p3<-ggplot(all,aes(x=risk_score,y=Age,col=AgeGroup))+
geom_point(na.rm = TRUE)+
stat_smooth(method =  lm)+
labs(x='Risk Score',y='Age')+
annotate("text", x = loc_x, y = loc_y, size = 4, colour = "black",
label=paste0("Adult r=",round(r2,3),"  p=",signif(p2,3),"\n",
"Child r=",round(r1,3),"  p=",signif(p1,3)
))+
theme_bw()+
theme(legend.position="right",
plot.title = element_text(hjust = 0.5))


ggsave("Risk_plt_SG10K_Age.pdf",p3,height=4,width=5)
ggsave("Risk_plt_SG10K_Age.png",p4,height=4,width=5)


p3_x<-ggplot(all,aes(x=risk_score,fill=AgeGroup))+
geom_density(alpha=0.8)+
theme_bw()+
theme(
axis.ticks=element_blank(),
text=element_blank(),
legend.position="none")

p3_y<-ggplot(all,aes(x=Age,fill=AgeGroup))+
geom_density(alpha=0.8)+
theme_bw()+
theme(
axis.ticks=element_blank(),
#text=element_blank(),
legend.position="none")+
coord_flip()

pdf("hist_risk_score.pdf",height=4,width=4)
hist(all$risk_score)
dev.off()

pdf("hist_age.pdf",height=4,width=4)
hist(all$Age)
dev.off()



st<-table(all$Group,all$Gender)
st
p.value<-fisher.test(st)$p.value
st[1,]<-st[1,]/sum(st[1,])
st[2,]<-st[2,]/sum(st[2,])
st1<-as.data.frame(st)
colnames(st1)<-c("Group","Gender","Freq")
write.csv(st,"Stat_Group_Gender.csv")

p4<-ggplot(st1,aes(x=Group,y=Freq,fill=Gender))+
geom_bar(stat="identity",position="stack")+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme_bw()+
xlab("Risk Group")+ylab("Percentage of samples")+
ggtitle(paste0("Fisher's P = ",round(p.value,3)))+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Freq,3)*100,"%")),
hjust = "center",color="white")+
theme(text=element_text(hjust = 0.5),
legend.position="right",
plot.title = element_text(hjust = 0.5))

ggsave("Risk_plt_SG10K_Gender.pdf",p4,height=4,width=4)
ggsave("Risk_plt_SG10K_Gender.png",p4,height=4,width=4)


plt3<-plot_grid(p3,p4,rel_widths=c(1.3,1))
plt3

ggsave("Risk_plt_SG10K_Gender_Age.pdf",plt3,height=4,width=8)
ggsave("Risk_plt_SG10K_Gender_Age.png",plt3,height=4,width=8)



plt4<-ggplot(all,aes(x=Group,y=risk_score,fill=Group))+
geom_violin(width=0.7)+
geom_boxplot(width=0.05,fill="white")+
facet_wrap(.~Ethnicity,scales="free_y")+
labs(x='Risk Group', y= "Risk score")+
geom_signif(comparisons = list(c("High Risk","Low Risk")),
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
theme_bw()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
legend.position='none',
axis.text.x=element_blank(),
plot.title = element_text(size=15,hjust = 0.5))

ggsave("Risk_plt_SG10K_Ethnicity.pdf",plt4,height=4,width=4)
ggsave("Risk_plt_SG10K_Ethnicity.png",plt4,height=4,width=4)