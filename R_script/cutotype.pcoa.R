infile = "~species.profile"
inmap ="~/Mappinp_file.txt"
incluster='./forehead_cluster.txt'
filename='./FH.Bray.PcoA2.pdf'
ino = "cluster"
inm = "pcoa"
location='forehead'
inc = 2

suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(clusterSim))
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggpubr))

#judge location
if(location=='forehead'){
  n1='FhSampleID'
  n2='FH'
}
if(location=='cheek'){
  n1='CkSampleID'
  n2='CK'
}
if(location=='nose'){
  n1='NsSampleID'
  n2='NS'
}

#read file
cluster<-read.table(incluster,header = T,check.names = F,sep='\t',stringsAsFactors = F)

data=read.csv(infile, header=T, check.names = F, sep="\t",row.names=1,stringsAsFactors = F)

map<-read.csv(inmap,header=T, check.names = F, sep="\t",stringsAsFactors = F)
mapping<-read.csv(inmap,header=T, check.names = F, sep="\t",stringsAsFactors = F)

data<-data[,colnames(data)%in%map$Sample]
map<-map[map$Sample%in%colnames(data),]
sum_type2<-table(map$type2)
sum_type2<-as.data.frame(sum_type2)
sum_type2<-sum_type2[sum_type2$Freq==3,]
map<-map[map$type2%in%sum_type2$Var1,]

map<-data.frame(SampleID=map[,1],Group=map[,2],type2=map[,3])

map<-map[map$Group==location,]
map<-map[map$SampleID %in% colnames(data),]
data<-data[,colnames(data)%in%map$SampleID]
data<-data[rowSums(data)!=0,]


#####################Bray-Curtis###########
library(vegan)
dist.bray<-vegdist(t(data))
pcoa <- cmdscale(dist.bray, k=10,eig = TRUE)


pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)


sample_site <- data.frame({pcoa$point})[1:2]
sample_site$names <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')

#
rownames(cluster)<-cluster[,n1]
cluster<-cluster[colnames(data),]
group<-data.frame(ID=colnames(data),cluster=cluster[,n2])

#add cluster information
sample_site<-sample_site[group$ID,]
sample_site$cluster<-group$cluster

adonis1<-adonis(dist.bray~ sample_site$cluster,permutations = 999,method = "bray")

sample_site$cluster<-factor(sample_site$cluster)
library(ggplot2)

#pcoa_plot
pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2, group = cluster,color=cluster)) +
  theme_classic()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  geom_point(size = 1.5) + #change point size
  scale_color_manual(values = brewer.pal(9,"Set1"),
                     breaks=c("1","2"),
                     labels=c("M.osloensis","C.acnes")) + #change color
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title=element_blank(),
        legend.position = c('none')
  )+
  labs(x = paste('PCoA1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA2: ', round(100 * pcoa_eig[2], 2), '%')) 

d<-ggplot(sample_site,aes(PCoA1,group=cluster,color=cluster,fill=cluster))+
  geom_bar()+
  geom_density(alpha=0.6)+
  scale_color_manual(values=brewer.pal(3,'Set1'))+
  scale_fill_manual(values=brewer.pal(3,'Set1'))+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ),
        legend.position = c('none'))+
  xlab('')+
  ylab('')

b<-ggplot(sample_site,aes(PCoA2,group=cluster,color=cluster,fill=cluster))+
  geom_bar()+
  geom_density(alpha=0.6)+
  scale_color_manual(values=brewer.pal(3,'Set1'))+
  scale_fill_manual(values=brewer.pal(3,'Set1'))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ),
        legend.position = c('none'))+
  xlab('')+
  ylab('')+coord_flip()

dist.bray<-as.matrix(dist.bray)
for(i in 1:nrow(dist.bray)){
  for(j in 1:ncol(dist.bray)){
    if(j>=i){
      dist.bray[i,j]=""
    }
  }
}
dist.bray<-reshape2::melt(dist.bray)
dist.bray<-dist.bray[dist.bray$value!="",]
dist.bray$value<-as.double(as.character(dist.bray$value))
dist.bray<-merge(dist.bray,group,by.x='Var1',by.y='ID',all.x=T)
dist.bray<-merge(dist.bray,group,by.x='Var2',by.y='ID',all.x=T)
a<-mean(dist.bray[dist.bray$cluster.x!=dist.bray$cluster.y,]$value)

dist.bray<-dist.bray[dist.bray$cluster.x==dist.bray$cluster.y,]
d2<-aggregate(dist.bray$value,by=list(cluster1=dist.bray$cluster.x,cluster2=dist.bray$cluster.y),mean)
d3<-aggregate(dist.bray$value,by=list(cluster1=dist.bray$cluster.x,cluster2=dist.bray$cluster.y),quantile)
d3<-d3$x[,2:4]
d2<-cbind(d2,d3)
m<-ggplot(d2,aes(x=factor(cluster1),y=x,group=factor(cluster1)))+
  geom_bar(stat='identity',aes(fill=factor(cluster1)),width = .7)+
  scale_fill_manual(values = brewer.pal(9,"Set1"),
                    breaks=c("1","2"),
                    labels=c("M.osloensis","C.acnes")) +
  geom_errorbar(aes(ymin=d2$`25%`, ymax=d2$`75%`), width=.1)+
  theme_classic()+
  geom_hline(aes(yintercept=a), colour="#990000")+
  guides(fill=F)+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )+
  xlab("")+ylab("")



p <- ggplotGrob(pcoa_plot)
d <- ggplotGrob(d)
b <- ggplotGrob(b)
m <- ggplotGrob(m)

pdf(filename,width = as.double(6), height = as.double(5))
grid.arrange(d,m,p,b,ncol=2,nrow=2,widths=c(4,1),heights = c(1,4))

dev.off()

