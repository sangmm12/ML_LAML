rm(list=ls())
file_dir = "~/diff/ml-laml"
setwd(dir = file_dir)
library(ggplot2)
library(ggpubr)
library(stringr)
need_list = c('rfrisk')

for(file_name in need_list)
{

  #tumor_list = c('ACC','BLCA','BRCA','CESC','COAD', 'CHOL', 'ESCA', 'GBM', 'HNSC', 'KIRP', 'KIRC', 'KICH', 'LGG', 'LUAD', 'LUSC', 'LIHC', 'LAML', 'OV', 'PRAD', 'PAAD', 'PCPG', 'READ', 'STAD', 'SKCM', 'THCA', 'TGCT', 'UCEC', 'UCS')
  
  tumor_list = c('LUAD')
  
  final_tumor_list = c()
  
  my_data = read.csv('~/diff/mac-m6a/GSE72094exp.csv',header=T,check.names=F)
  
  other_data = read.csv(paste(file_name,".csv",sep=''),header=T,check.names=F)
  
  p_value_csv <- data.frame()
  
  for(name in tumor_list)
  {
  
    dat <- data.frame(check.names = F)
    
    other_file <- other_data[grep(name,other_data$CODE),]
    
    if (length(rownames(other_file))==0)
    {
        next
    }
    
    other_file <- other_file[!duplicated(other_file$SampleName),]
    
    rownames(other_file) = gsub('\n','',other_file$SampleName)
    
    other_file <- subset(other_file,select=-c(SampleName,CODE))
    
    other_file <- other_file[grep("Tumor",other_file$Group),]
    other_file <- subset(other_file,select=-c(Group))
    
    exp_file <- my_data[grep(name,my_data$CODE),]
    
    exp_file <- exp_file[grep("Tumor",exp_file$Group),]
    
    exp_file <- exp_file[!duplicated(exp_file$SampleName),]
    
    rownames(exp_file) = gsub('\n','',exp_file$SampleName)
    
    exp_file <- subset(exp_file,select=-c(Group,CODE,SampleName))
    
    gene_list <- colnames(exp_file)
    
    all_name <- names(which(table(c(rownames(other_file),rownames(exp_file)))==2))
    
    if (length(all_name)==0)
    {
        next
    }
    
    for(gene_name in gene_list)
    {
        for(i in all_name)
        {
            dat <- rbind(dat,c(gene_name,other_file[match(i,rownames(other_file)),],exp_file[c(gene_name)][match(i,rownames(exp_file)),]))
        }
    
    }
    
    colnames(dat) <- c("Gene","Group","value")
    
    dat[,3] = as.numeric(dat[,3])
    
    
    dat <- na.omit(dat)
    xx <-compare_means(value ~ Group, data = dat, group.by = "Gene",method = "anova")
    
    
    p_value <- as.matrix(xx$p)
    
    p_value[is.na(p_value)] <- 1
    
    
    
    final_tumor_list <- append(final_tumor_list,name)
    print(name)
    if(length(final_tumor_list)!=1){p_value_csv <- cbind(p_value_csv,as.data.frame(p_value))}else{p_value_csv <- as.data.frame(p_value)}
    
    
    pdf(paste(file_name,"-",name,".pdf",sep=''),width=length(unique(dat[,1]))-5,height = 8)
    
    p <- ggboxplot(dat, x = "Gene", y = "value",
                  color = "Group", palette = "jama",
                  add = c("jitter"),x.text.angle=60)
    

    p <- p + xlab("")+ylab("Gene Expression(log(x+1))")
    p <- p + theme(axis.text = element_text(size = 30),axis.title=element_text(size=30))
    print(p + stat_compare_means(aes(group = Group), label = "p.signif",method = "anova"))
    
    dev.off()
  }
  rownames(p_value_csv) <- gene_list
  colnames(p_value_csv) <- final_tumor_list
  write.csv(p_value_csv,file=paste("p.csv",sep=''),quote=F)
}

