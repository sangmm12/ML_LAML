setwd('~/machine_learning/cox')


library(dplyr)
library(survival)
library(data.table)
store_data = list()

for(sample_name in c('tcga','GSE12417','GSE37642','GSE146173','GSE106291'))
{

dat <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T,row.name=1,check.names = F)

all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)

output <- data.frame()

for(gene_name in colnames(dat))
{
  dat1 <- dat[which(colnames(dat)==gene_name)]
  
  dat1$SampleName = rownames(dat1)

  erged_df <- merge(dat1, all_sur_data, by = "SampleName")
  
  erged_df <- erged_df[,-c(1)]
  
  cox_model <- coxph(Surv(Time,Status) ~ as.numeric(erged_df[,1]),data=erged_df)
  
  lower = round(summary(cox_model)$conf.int[,3],2)
  upper = round(summary(cox_model)$conf.int[,4],2)
  hr = round(summary(cox_model)$conf.int[,1],2)
  p = summary(cox_model)$coefficients[,5]
  if(is.na(p)||p>0.05||hr<1)
  {
    next
  }
  ci = paste(hr,'(',lower,'|',upper,')',sep='')
  temp_data <- data.frame(genename = c(gene_name),pvalue = c(p),HR=c(hr),Lower=c(lower),Upper=c(upper),CI=c(ci))
  output <- rbind(output,temp_data)
  eval(parse(text = paste('store_data$',sample_name,"<-append(store_data$",sample_name,",'",gene_name,"')",sep='')))
}
}