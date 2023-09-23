
library(dplyr)
library(corrplot)

setwd("~/cor_tri")

myColor <- colorRampPalette(c("#B7E731", "white", "#BD2989"))(1000)

my_data <- read.csv('~/diff/data/CIBER.csv',header=T,check.names=F)
my_data <- my_data[grep('LAML',my_data$CODE),]
my_data <- subset(my_data,select=-c(Group,CODE,SampleName))

res <- cor(my_data)
#method = "square"
pdf("1.pdf",width=length(colnames(my_data))/3,height=length(colnames(my_data))/3)
corrplot(res, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 90,col=myColor)
dev.off()