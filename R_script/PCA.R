
suppressPackageStartupMessages(library(ade4))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(vegan))
suppressPackageStartupMessages(library(ggpubr))


AD<-read.csv('./species.txt',header=T,sep='\t',row.names = 1,stringsAsFactors = F,check.names = F,quote="")
HMP<-read.csv('./hmpspcies.txt',header=T,sep='\t',row.names = 1,stringsAsFactors = F,check.names = F,quote="")
map<-read.table('./mapping.txt',,header=T,sep='\t',stringsAsFactors = F,check.names = F)
map<-map[map$Type=='Ac',]
AD<-AD[rownames(AD)%in%rownames(HMP),]
HMP<-HMP[rownames(HMP)%in%rownames(AD),]
HMP<-HMP[rownames(AD),]
HMP<-cbind(HMP,AD)
HMP<-HMP[,colnames(HMP)%in%map$Sample]

dat<-HMP
data1<-data.frame('Sample'=map$Sample,'Group'=map$Project,'Condition'=map$Condition)
dat<-dat[,map$Sample]
rownames(dat)<-gsub('.*;','',rownames(dat))
rownames(data1)<-data1$Sample
data = dat
colnames(data1)[2] ="Group"
data1$Group<-factor(data1$Group)

HMP$mean<-rowSums(HMP)/ncol(HMP)
HMP<-HMP[order(HMP$mean,decreasing = T),]
rownames(HMP)<-gsub('.*;','',rownames(HMP))
top<-rownames(HMP)[1:12]

incw = 1.5
inlength = 1
inlift = "Axis1"
inright = "Axis2"
newpalette<-c('#4458BB','#D4789C',
              '#82B6FC','#AA9CD6','#F1CED5')


data<-t(sqrt(data))

data.dudi <- dudi.pca(data, center=TRUE, scale=F, scan=F, nf=10)
data2 <- data.dudi$li
incw = as.double(incw)
inlength = as.double(inlength)
classified_c = as.character(unique(data1[,2]))

#Microbiome
adonis1<-adonis(t(dat) ~ data1[,2],permutations = 999,method = "bray")


phenotype <- data1[,2]
f=classified_c
Type <- factor(phenotype,levels=f)
m = data.dudi$li
n = data.dudi$c1

lift_data = m[as.character(inlift)]
row.names(lift_data) = gsub(pattern = "[.]",replacement = "-",row.names(lift_data))
right_data = m[as.character(inright)]
row.names(right_data) = gsub(pattern = "[.]",replacement = "-",row.names(right_data))
data.dudi$li = cbind(lift_data,right_data)
num1 = substring(as.character(inlift),5,6)
num2 = substring(as.character(inright),5,6)
num1_data = n[paste("CS",num1,sep = '')]
num2_data = n[paste("CS",num2,sep = '')]
data.dudi$c1 = cbind(num1_data,num2_data)

right_data_box= cbind(data1[,2],right_data)
colnames(right_data_box)[1] = "Group" 
lift_data_box = cbind(data1[,2],lift_data)
colnames(lift_data_box)[1] = "Group" 

#png(paste(infile,"PC",num1,"-","PC",num2,".png",sep = "_"), width = 1000, height = 768, res = 100)
#png("2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.Disease.PCA.png", width = 768, height = 768, res = 72)
x1 <- min(data.dudi$li[,1]) - 0.3
y1 <- min(data.dudi$li[,2]) - 0.3
x2 <- max(data.dudi$li[,1]) + 0.3
y2 <- max(data.dudi$li[,2]) + 0.3
bb <- head(data.dudi$c1[order(sqrt((data.dudi$c1[,1])^2+(data.dudi$c1[,2])^2),decreasing=T),],n=7L)[1:7,]
#bb<-data.dudi$c1[rownames(data.dudi$c1)%in%top,]
rownames(bb) <- gsub("^X", "", rownames(bb))
rownames(bb) <- gsub("\\S+o__", "o__", rownames(bb))
cutoff <- (x2-0.3) / abs(bb[1,1]) * inlength
d2 <- data.frame(X=bb[1:dim(bb)[1],1]*cutoff, Y=bb[1:dim(bb)[1],2]*cutoff, LAB=rownames(bb)[1:dim(bb)[1]])
#d2[[3]]<- gsub('.*f__','f__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*o__','o__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*c__','c__',as.character(d2[[3]]))
d2[[3]] <- gsub('.*p__','p__',as.character(d2[[3]]))

d <- data.dudi$li
eig <- ( data.dudi$eig / sum( data.dudi$eig ) )

#track <- read.table("2.beta_div/plots/pca/2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.track", head=T, sep="	")
#points <- read.table("2.beta_div/plots/pca/2017.Jun11.Throat.Stool.16S.otu.Tonsil.prof.points", head=T, sep="	")

ggdata <- data.frame(d)

data1<-data1[rownames(d),]

p<-ggplot(ggdata) +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(x=d[,1], y=d[,2], color=Type,shape=factor(data1$Condition,levels=c('Control','Case'))), size=3) +
  geom_text_repel(data=d2, aes(x=X, y=Y, label=LAB),
                  family="Helvetica", fontface="italic", size=3) +
  geom_segment(data=d2, aes(x=0, y=0, xend=X, yend=Y),
               arrow = arrow(length = unit(0.3, "cm")), size=0.8, alpha=0.5)+
  scale_color_manual(values=newpalette)+
  scale_fill_manual(values=newpalette)+
  guides(color=guide_legend(colnames(data1)[2]),
         fill=guide_legend(colnames(data1)[2]),
         shape=guide_legend(colnames(data1)[2]) ) +
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position=c(0.9,0.8),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=9
        ),
        legend.title=element_text(
          family="Helvetica",
          colour="black",
          size=6
        ),
        legend.text=element_text(
          family="Helvetica",
          colour="black",
          size=10
        )
  )+xlim(min(d[,1], 0)*incw, max(d[,1])*incw)+ylim(min(d[,2], 0)*incw, max(d[,2])*incw)

p <- ggplotGrob(p)


right_data_box$Group = factor(right_data_box$Group, levels=f)
d <- ggplot(right_data_box,aes(x = Group,y = right_data_box[,2],fill = Group))+
  geom_boxplot(width = 0.5)+
  theme_bw()+theme(panel.grid =element_blank())+
  scale_fill_manual(values=newpalette,breaks =f)+
  stat_compare_means()+
  guides(fill=FALSE)+theme(axis.text.x = element_blank())+
  theme(axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  ylim(min(right_data_box[,2], 0)*incw, max(right_data_box[,2])*incw)+
  xlab("")+ylab(paste("PC",num2," (",round(eig[as.numeric(num2)]*100,2),"%)",sep=""))
lift_data_box$Group = factor(lift_data_box$Group, levels=f)
b<- ggplot(lift_data_box,aes(x = Group,y = lift_data_box[,2],fill = Group))+
  geom_boxplot(width = 0.5)+
  theme_bw()+
  theme(panel.grid =element_blank())+coord_flip()+
  guides(fill=FALSE)+
  theme(axis.text.y = element_blank())+
  theme(axis.ticks = element_blank(),
        text=element_text(
          family="Helvetica",
          face="bold",
          colour="black",
          size=12
        ))+
  scale_fill_manual(values=newpalette)+
  stat_compare_means()+
  ylim(min(lift_data_box[,2], 0)*incw, max(lift_data_box[,2])*incw)+
  xlab("")+ylab(paste("PC",num1," (",round(eig[as.numeric(num1)]*100,2),"%)",sep=""))
a<-ggplot()+theme_bw()+theme(panel.border = element_blank(),panel.grid =element_blank(),
                             axis.text = element_blank(),axis.title = element_blank(),
                             axis.ticks = element_blank())+
  annotate("text", x=1, y=40, label=paste("P.value =",round(adonis1[[1]][6][[1]][1],4),'\n',
                                          "R2      =",round(adonis1[[1]][5][[1]][1],4)), size=3.5)

a <- ggplotGrob(a)
d <- ggplotGrob(d)
b <- ggplotGrob(b)


pdf('./DiffCountry.PCA.pdf',width = as.double(9), height = as.double(6))
grid.arrange(d,p,a,b,ncol=2,widths=c(1,4),heights = c(4,1))

dev.off()

