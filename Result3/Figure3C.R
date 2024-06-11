setwd("/p300s/jiapl_group/zhuxl/Microbiome/Cor")


library(ggplot2)
library(dplyr)
library(ggcor)
library(linkET)


data<-read.table("mpa_table.xls",header=T,sep="\t",row.names=1)
data[1:4,1:4]
colnames(data)<-gsub("X","",colnames(data))
species<-t(data)
species[1:4,1:4]


mat<-read.csv("mat_manhattan.csv",row.names=1)
mat[1:4,1:4]
colnames(mat)<-gsub("X","",colnames(mat))
group<-read.table("group.txt",header=T,sep="\t")
group<-group[order(group$Group),]
mat1<-mat[,group$SampleID]
table(colnames(mat1)==group$SampleID)
mat2<-mat1[rownames(species),]
table(rownames(mat2)==rownames(species))

write.table(mat2,"mat2_order.tsv",sep="\t",quote=F)
write.table(species,"species_order.tsv",sep="\t",quote=F)


mantel <- mantel_test(mat2,species,
spec.select = list(Group1 = 1:30,Group2 = 31:373,
Group3 = 374:964,Group4 = 965:993)) %>% 
mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

save(mantel,file="mantel.RData")

mantel1<-mantel[mantel$p.value<0.05,]


plt<-quickcor(species, type = "upper") +
geom_square() +
anno_link(aes(colour = pd, size = rd), data = mantel1) +
scale_size_manual(values = c(0.5, 1, 2)) +
scale_fill_gradientn(colors = c('#053061', '#68A8CF', 'white', '#F7B394', '#67001F'), limits = c(-1, 1)) +
scale_colour_manual(values = c('#56B4E9', '#E69F00', '#999999')) +
guides(size = guide_legend(title = "Mantel's r",override.aes = list(colour = "grey35"), order = 2),
colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
fill = guide_colorbar(title = "Pearson's r", order = 3))+
#add_diag_label(angle = 45) +  #如果想在对角线显示环境变量标签
#remove_axis('all') +
theme(legend.key = element_blank())

ggsave("Group_Species_Cor.pdf",plt,height=10,width=10)