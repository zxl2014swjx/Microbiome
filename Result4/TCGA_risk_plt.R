setwd("D:/Microbiome/5.TCGA")

taxa<-read.table("TCGA_microbiome.xls",,header=T,sep="\t",row.names=1)
dim(taxa)
taxa<-taxa[,-1]
taxa[1:4,1:4]
taxa[is.na(taxa)]<-0
taxa<-taxa[,colSums(taxa)!=0]
write.csv(taxa,"taxa_input.csv")


H<-read.csv("Genus_Group_sum.csv")
head(H)
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


#G3_shannon<-vegan::diversity(t(G3_metagenome),'shannon')
#G2_shannon<-vegan::diversity(t(G2_metagenome),'shannon')


G3_shannon<-microbiome::alpha(G3_metagenome,index=c("diversity_shannon"))
G2_shannon<-microbiome::alpha(G2_metagenome,index=c("diversity_shannon"))


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
colnames(mer)<-c("G3_shannon","G2_shannon","R_G3","R_G2",
"psi_G3","psi_G2","risk_score","sample")
head(mer)
write.csv(mer,"TCGA_Risk_Score.csv",row.names=F)

###############


setwd("D:/Microbiome/5.TCGA")

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
BCI<-BCI+10
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
data<-read.xlsx("TCGA_Risk_Score.xlsx",sheet=1)
dim(data)
head(data)

all<-data

table(all$Group)
summary(all$risk_score)


score<-all$risk_score
pdf("Density_Risk_Group.pdf",height=4,width=4)
r <- hist(score, col = "lightblue", border = "grey50",
freq = F)
text(r$mids, r$density, r$counts, pos=3,col = "black")
lines(density(score), col = "red")
dev.off()

all[all$risk_score>0,]$Group<-"High Risk"
all[all$risk_score<=0,]$Group<-"Low Risk"
table(all$Group)

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

ggsave("Risk_plt_TCGA.pdf",plt,height=4,width=8)
ggsave("Risk_plt_TCGA.png",plt,height=4,width=8)


##########


need<-all[,c(1,2,9:30)]
need1<-melt(need,id.vars=1:4)
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

ggsave("Risk_plt_TCGA_diversity.pdf",plt2,height=12,width=15)
ggsave("Risk_plt_TCGA_diversity.png",plt2,height=12,width=15)


################
st<-compare_means(risk_score~Group, data=all,group.by= "Cancer")


plt3<-ggplot(all,aes(x=Cancer,y=risk_score,fill=Group))+
geom_boxplot(width=0.4)+
labs(x='Risk Group', y= "Risk score")+
geom_signif(comparisons = list(c("High Risk","Low Risk")),
step_increase = 0.1,map_signif_level = F,test = wilcox.test)+
theme_bw()+
scale_fill_manual(values = c("#CB3425","#3F5688"))+
theme(legend.title=element_blank(),
legend.position='top',
text=element_text(size=15,hjust = 0.5),
axis.text.x=element_text(angle=90),
plot.title = element_text(size=15,hjust = 0.5))


ggsave("Risk_plt_TCGA_Group.pdf",plt3,height=5,width=7)
ggsave("Risk_plt_TCGA_Group.png",plt3,height=5,width=7)



