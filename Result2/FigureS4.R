setwd("D:/Microbiome/Disease_micorbiome_pairs")

all<-read.table("Phenotype_microbiome_relationship.xls",header=T,sep="\t")
colnames(all)[21:38]


table(all$T2D)
table(all$T2D.drug)
table(all$T2D.family.history)
table(all$hypertension)
table(all$hypertension.drug)
table(all$hypertension.family.history)
table(all$high.blood.lipid)
table(all$high.blood.lipid.drug)
table(all$high.blood.lipid.family.history)
table(all$heart.disease)
table(all$heart.disease.drug)
table(all$heart.disease.family.history)
table(all$cerebral.hemorrhag)
table(all$cerebral.hemorrhag.drug)
table(all$cerebral.hemorrhag.family.history)
table(all$cerebral.thrombosis)
table(all$cerebral.thrombosis.drug)
table(all$cerebral.thrombosis.family.history)

#####################

stat1<-table(all$T2D,all$T2D.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("T2D","T2D.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=T2D,y=Num,fill=T2D.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=T2D,y=Num,fill=T2D.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$T2D,all$T2D.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("T2D","T2D.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=T2D,y=Num,fill=T2D.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=T2D,y=Num,fill=T2D.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}

if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_T2D.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_T2D.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_T2D.pdf",plt1,height=5,width=4)
}
rm(p1,p2)



#####################

stat1<-table(all$hypertension,all$hypertension.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("hypertension","hypertension.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=hypertension,y=Num,fill=hypertension.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=hypertension,y=Num,fill=hypertension.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$hypertension,all$hypertension.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("hypertension","hypertension.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=hypertension,y=Num,fill=hypertension.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=hypertension,y=Num,fill=hypertension.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_hypertension.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_hypertension.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_hypertension.pdf",plt1,height=5,width=4)
}
rm(p1,p2)

#####################

stat1<-table(all$high.blood.lipid,all$high.blood.lipid.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("high.blood.lipid","high.blood.lipid.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=high.blood.lipid,y=Num,fill=high.blood.lipid.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=high.blood.lipid,y=Num,fill=high.blood.lipid.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$high.blood.lipid,all$high.blood.lipid.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("high.blood.lipid","high.blood.lipid.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=high.blood.lipid,y=Num,fill=high.blood.lipid.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=high.blood.lipid,y=Num,fill=high.blood.lipid.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_high.blood.lipid.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_high.blood.lipid.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_high.blood.lipid.pdf",plt1,height=5,width=4)
}
rm(p1,p2)

################

stat1<-table(all$heart.disease,all$heart.disease.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("heart.disease","heart.disease.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=heart.disease,y=Num,fill=heart.disease.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=heart.disease,y=Num,fill=heart.disease.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$heart.disease,all$heart.disease.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("heart.disease","heart.disease.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=heart.disease,y=Num,fill=heart.disease.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=heart.disease,y=Num,fill=heart.disease.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_heart.disease.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_heart.disease.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_heart.disease.pdf",plt1,height=5,width=4)
}
rm(p1,p2)

###################
stat1<-table(all$cerebral.hemorrhag,all$cerebral.hemorrhag.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("cerebral.hemorrhag","cerebral.hemorrhag.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=cerebral.hemorrhag,y=Num,fill=cerebral.hemorrhag.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=cerebral.hemorrhag,y=Num,fill=cerebral.hemorrhag.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$cerebral.hemorrhag,all$cerebral.hemorrhag.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("cerebral.hemorrhag","cerebral.hemorrhag.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=cerebral.hemorrhag,y=Num,fill=cerebral.hemorrhag.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("Drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=cerebral.hemorrhag,y=Num,fill=cerebral.hemorrhag.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_cerebral.hemorrhag.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_cerebral.hemorrhag.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_cerebral.hemorrhag.pdf",plt1,height=5,width=4)
}
rm(p1,p2)
###########################

stat1<-table(all$cerebral.thrombosis,all$cerebral.thrombosis.drug)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("cerebral.thrombosis","cerebral.thrombosis.drug","Num")
if(pval1<0.01){
p1<-ggplot(stat1,aes(x=cerebral.thrombosis,y=Num,fill=cerebral.thrombosis.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p1<-ggplot(stat1,aes(x=cerebral.thrombosis,y=Num,fill=cerebral.thrombosis.drug))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("drug")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


stat1<-table(all$cerebral.thrombosis,all$cerebral.thrombosis.family.history)
pval1<-fisher.test(stat1)$p.value
for(i in 1:dim(stat1)[1]){
stat1[i,]<-stat1[i,]/sum(stat1[i,])
}
stat1<-stat1[order(stat1[,2]),]
label<-rownames(stat1)
stat1<-as.data.frame(stat1)
colnames(stat1)<-c("cerebral.thrombosis","cerebral.thrombosis.family.history","Num")
if(pval1<0.01){
p2<-ggplot(stat1,aes(x=cerebral.thrombosis,y=Num,fill=cerebral.thrombosis.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",format(pval1,scientific=TRUE,digit=3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}
if(pval1>0.01){
p2<-ggplot(stat1,aes(x=cerebral.thrombosis,y=Num,fill=cerebral.thrombosis.family.history))+
geom_bar(stat="identity",position="fill",width=0.5)+
theme_bw()+
scale_fill_nejm()+
ggtitle(paste0("Fisher test p = ",round(pval1,3)))+
theme(legend.position="none",
text=element_text(size=15),
axis.text.x=element_text(angle=0,vjust=0.5),
plot.title = element_text(hjust = 0.5))+
ylab("family history")+
scale_y_continuous(expand=c(0,0))+
geom_text(position=position_stack(0.5),
aes(label = paste0(round(Num,3)*100,"%")),
hjust = "center",size=3,color="black")+
scale_x_discrete(limit=label)
}


if(exists("p1") & exists("p1")){
plt1<-plot_grid(p1,p2,nrow=1)
ggsave("Bar_cerebral.thrombosis.pdf",plt1,height=5,width=7)
}
if(exists("p1") & !exists("p2")){
plt1<-plot_grid(p1,nrow=1)
ggsave("Bar_cerebral.thrombosis.pdf",plt1,height=5,width=4)
}
if(exists("p2") & !exists("p1")){
plt1<-plot_grid(p2,nrow=1)
ggsave("Bar_cerebral.thrombosis.pdf",plt1,height=5,width=4)
}
rm(p1,p2)