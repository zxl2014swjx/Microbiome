# Figure S1B


data<-read.table("Bam_Q30.xls",header=T,sep="\t")


p1<-ggplot(data,aes(x=Q30))+
geom_density(alpha=0.8)+
theme_bw()+
theme(legend.title=element_blank(),
text=element_text(size=15,hjust = 0.5),
legend.position="top")


ggsave("Kraken_Filtered_Stat.png",p1,height=4,width=4)
ggsave("Kraken_Filtered_Stat.pdf",p1,height=4,width=4)


pdf("Q30_density.pdf",height=8,width=5)
r <- hist(Q30, col = "lightblue", border = "grey50",
freq = F, breaks = 10,ylim=c(0,1))
text(r$mids, r$density, r$counts, pos=3,col = "black")
lines(density(Q30), col = "red")
dev.off()