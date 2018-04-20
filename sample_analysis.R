require(ggplot2)
require(plyr)
require(grid)

source("marpol.R")
source("polarCoord.R")
source("angular-analysis.R")

# Load an ERES file
eres<-read.csv('yipf_golgi1_golgi-example.csv')

# Load Nuclei file
nuclei<-read.csv('yipf_golgi1_Nuclei_example.csv')

# Join both datasets to add nuclei coordinates to the table
joined<-merge(eres,nuclei,by.x=c("ImageNumber","Parent_Cells","Metadata_Date"),by.y=c("ImageNumber","ObjectNumber","Metadata_Date"))

# Create new dataframe with coordinates only
data<-data.frame(x=joined$Location_Center_X.x-joined$Location_Center_X.y, 
                 y=joined$Location_Center_Y.x-joined$Location_Center_Y.y, 
                 cell=paste(joined$ImageNumber,joined$Parent_Cells,sep=""),
                 trt=joined$Metadata_Date)

data.polar<-data.frame(cart2pol(data$x,data$y),cell=data$cell,trt=data$trt)

# normalise angles by subtracting average angle
normAngs<-ddply(data.polar,.(trt,cell),transform,n=(angle-meanAngle(angle)) %% (2*pi),.progress="text")

# histograms function

# angular histograms
angHisto<-ddply(normAngs,.(trt,cell),function(x) calcHist1(x$n),.progress="text")

# calculate histograms
histoDif<-ddply(angHisto,.(trt,cell),function(x) calcDiff(x[-1:-3]),.progress="text")

# calculate average histoDif
means<-ddply(histoDif,.(trt),summarise,mean=mean(V1),sd=sd(V1),se=sd(V1)/sqrt(length(V1)),N=length(V1))
write.csv(means,"means.csv")

# Do a statistical test. Most of the groups have normal distribution, but some don't. Therefore we perform Wilcoxon test with FDR
# boxplot shows some skewness in several samples
#ggplot(histoDif, aes(trt,V1))+geom_boxplot()

# do a wilcox test and select only comparisons to Neg
pvals<-pairwise.wilcox.test(histoDif$V1,histoDif$trt,p.adjust.method = "BH")
pvals<-data.frame(rows=rownames(as.data.frame(pvals$p.value)),stack(as.data.frame(pvals$p.value)))
pvals<-na.omit(pvals)
pvals.p<-pvals[pvals$rows=="neg" | pvals$ind=="neg",]
write.csv(pvals.p,"pvalues.csv")

# Plot the polar distribution scores
# get star offset and create labels per group
# Symbol  Meaning
# ns  P > 0.05
# *  P ≤ 0.05
# **	P ≤ 0.01
# ***	P ≤ 0.001
# ****	 P ≤ 0.0001

means$star<-c("","","","**","****","*","","*")


text_pt<-9
text_mm<-9*0.35
offset<-0.05

plot_polar<-ggplot(means,aes(x=trt,y=mean))+geom_bar(stat="identity",fill='grey50', width=0.5,size=0.1, colour="black")+
  geom_errorbar(aes(ymax=mean+se, ymin=mean-se), position='dodge',colour='black', width=0.1,size=0.1) +
  geom_text(aes(y=mean+offset,label=star),size=text_mm)
plot_polar + theme_bw() + 
  scale_y_continuous(limits=c(0,1.4),expand = c(0,0)) + labs(title="Polar distribution of Golgi fragments",x="siRNA Target", y="Polar distribution score") +
  theme(text=element_text(size=text_pt), plot.margin=unit(c(1,1,0,0),"mm"),plot.title=element_blank(),
        axis.text.x=element_text(size=text_pt, angle=40, hjust = 1)) +
  scale_x_discrete(labels=c("Neg","YIPF1","YIPF2","YIPF3","YIPF4","YIPF5","YIPF6","YIPF7")) +
  coord_cartesian(ylim=c(0.7,1))

ggsave("golgi-polar-paper_1.pdf",width=60,height=50,units="mm")