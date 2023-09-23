setwd("~/correlate")
rm(list = ls())
library(dplyr)
library(patchwork)
library(ggplotify)


tumor_list = c('LAML')

for (tumor_id in tumor_list)
{
  print(tumor_id)
  
  dat_input <- read.csv("data/machine_learning/hallmark.csv",check.names = F,header=T)
  
  dat_input <- dat_input[grep(tumor_id,dat_input$CODE),]
  
  dat_input <- dat_input[grep("Tumor",dat_input$Group),]
  
  dat_input <- subset(dat_input,select=-c(Group,CODE))
  
  dat_input <- dat_input[!duplicated(dat_input$SampleName),]
  
  rownames(dat_input) <- make.unique(sub('\n','',dat_input$SampleName))
  
  dat_input <- subset(dat_input,select=-c(SampleName))
  
  colscluster = length(colnames(dat_input))/3
  
  
  data_all_tumor <- read.csv(paste("data/machine_learning/gene_exp.csv",sep=''),check.names = F)
  
  dat_one_tumor <- data_all_tumor[grep(tumor_id,data_all_tumor$CODE),]
  
  if (length(colnames(dat_one_tumor))==0)
  {
    next
  }
  
  rownames(dat_one_tumor) <- make.unique(sub('\n','',dat_one_tumor$SampleName))
  
  dat_one_tumor <- dat_one_tumor[,!duplicated(colnames(dat_one_tumor))]
  
  dat_one_tumor <- dat_one_tumor[grep("Tumor",dat_one_tumor$Group),]
  
  dat_one_tumor <- subset(dat_one_tumor,select=-c(Group,SampleName,CODE))
  
  one_tumor_sample <- unlist(rownames(dat_one_tumor))
  
  all_name <- names(which(table(c(rownames(dat_input),one_tumor_sample))==2))
  
  dat_gene <- dat_one_tumor[match(all_name,rownames(dat_one_tumor)),]
  
  dat_im <- dat_input[match(all_name,rownames(dat_input)),]
  
  library(psych)
  data.corr <- corr.test(dat_gene, dat_im, method="pearson", adjust="fdr")
  
  data.r <- data.corr$r
  data.p <- data.corr$p
  
  library(pheatmap)
  
  
  getSig <- function(dc) {
    sc <- ' '
    if (dc < 0.0001) {sc <- '****'}
    else if (dc < 0.001){sc <- '***'}
    else if (dc < 0.01){sc <- '**'}
    else if (dc < 0.05) {sc <- '*'}
    else{sc <- ''
    }
    return(sc)
  }
  
  sig.mat <- matrix(sapply(data.p, getSig), nrow=nrow(data.p))
  str(sig.mat)
  
  
  paletteLength <- 1000
  # myColor <- colorRampPalette(c("#204E8C", "white", "#FFA61C"))(paletteLength)
  # myColor <- colorRampPalette(c("#ff8601", "white", "#044d69"))(paletteLength)
  # #gx
  # myColor <- colorRampPalette(c("#F59234", "white", "#20687D"))(paletteLength)
  # #cell
  # myColor <- colorRampPalette(c("#B7E731", "white", "#BD2989"))(paletteLength)
  # #mirna
  # myColor <- colorRampPalette(c("darkblue", "white", "red"))(paletteLength)
  # #hallmark
  myColor <- colorRampPalette(c("#00B060", "white", "#FF9000"))(paletteLength)
  # #KEGG
  #myColor <- colorRampPalette(c("#80E800", "white", "#1240AB"))(paletteLength)
  # #GO
  #myColor <- colorRampPalette(c("#E60042", "white", "#FFCD00"))(paletteLength)
  
  
  test <- data.r
  myBreaks <- seq(-1,1,length.out=paletteLength)
  
  
  pdf(paste(tumor_id,".pdf",sep=''),width = 30,height = 40)
  
  pheatmap(data.r, 
           color=myColor,
           breaks=myBreaks,
           clustering_method="average",border_color='grey',fontsize_row = 15,fontsize_col = 15,fontsize_number = 10, cluster_rows=F,cluster_cols=F,cellwidth = 25,cellheight = 25, display_numbers=sig.mat,fontsize=length(colnames(dat_input)))
  dev.off()
  print(rownames(data.r)[which(data.r == min(data.r), arr.ind = TRUE)[, 1]])
  print(colnames(data.r)[which(data.r == min(data.r), arr.ind = TRUE)[, 2]])
}
