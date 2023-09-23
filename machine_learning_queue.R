
setwd('~/machine_learning/new diagram')


library(dplyr)
library(tibble)
library(survival)
library(plsRcox)
library(superpc)
library(gbm)
library(CoxBoost)
library(survivalsvm)
library(BART)

coxphitermax = 3
SPCnfold = 3
final_result <- data.frame()

for(sample_name in c('tcga','GSE12417','GSE37642','GSE146173','GSE106291'))
{
  
  #####prepare for dataset
  result <- data.frame()
  all_sur_data <- read.table(paste('~/machine_learning/data/',sample_name,'/LAML_patient.txt',sep=''),header=T)
  gene_exp <- read.csv(paste('~/machine_learning/data/',sample_name,'/survival_out.csv',sep=''),header=T)
  pre_var <- colnames(gene_exp)
  pre_var <- pre_var[-1]
  mixed <- merge(gene_exp, all_sur_data, by = "SampleName")
  colnames(mixed) <- sub('Time','OS.time',colnames(mixed))
  colnames(mixed) <- sub('Status','OS',colnames(mixed))
  temp_row <- mixed$SampleName
  mixed<- mixed %>% 
    subset(select = -c(SampleName)) %>% 
    lapply(as.numeric) %>% 
    as.data.frame()
  mixed$SampleName <- temp_row
  mixed <- mixed %>% 
    select( SampleName,OS, OS.time, everything())
  seed <- 123456
  mixed <- na.omit(mixed)
  
  ### CoxBoost ####
  set.seed(seed)
  pen <- optimCoxBoostPenalty(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                              trace=TRUE,start.penalty=500,parallel = T)
  cv.res <- cv.CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                        maxstepno=500,K=10,type="verweij",penalty=pen$penalty)
  fit <- CoxBoost(mixed$OS.time,mixed$OS,as.matrix(mixed[,-c(1,2,3)]),
                  stepno=cv.res$optimal.step,penalty=pen$penalty)
  RS_COXBOOST <- data.frame(RS = as.numeric(predict(fit,newdata=mixed[,-c(1,2,3)],newtime=mixed[,3], newstatus=mixed[,2], type="lp")),name='CoxBoost')
  rs <- cbind(mixed[,2:3],RS = as.numeric(predict(fit,newdata=mixed[,-c(1,2,3)],newtime=mixed[,3], newstatus=mixed[,2], type="lp")))
  cc <- t(data.frame('CoxBoost'=summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### Survival SVM ####
  fit <- survivalsvm(Surv(OS.time,OS)~., data= mixed[,-c(1)], gamma.mu = 1)
  RS_SVM <- data.frame(RS=as.numeric(predict(fit, mixed[,-c(1)])$predicted),name='Survival SVM')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit, mixed[,-c(1)])$predicted))
  as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1])
  cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('Survival SVM')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### GBDT ####
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
             n.trees = 10000,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 6)
  best <- which.min(fit$cv.error)
  set.seed(seed)
  fit <- gbm(formula = Surv(OS.time,OS)~.,data = mixed[,-c(1)],distribution = 'coxph',
             n.trees = best,
             interaction.depth = 3,
             n.minobsinnode = 10,
             shrinkage = 0.001,
             cv.folds = 10,n.cores = 8)
  RS_GBDT <- data.frame(RS = as.numeric(predict(fit,mixed[,-c(1)],n.trees = best,type = 'link')),name='GBDT')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit,mixed[,-c(1)],n.trees = best,type = 'link')))
  cc <- data.frame(name = as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('GBDT')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### Supervised principal components ####
  data <- list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)])
  set.seed(seed)
  fit <- superpc.train(data = data,type = 'survival',s0.perc = 0.5) #default
  set.seed(seed)
  cv.fit <- superpc.cv(fit,data,n.threshold = 20,#default
                       n.fold = SPCnfold,
                       n.components=3,
                       min.features=5,
                       max.features=nrow(data$x),
                       compute.fullcv= TRUE,
                       compute.preval=TRUE)
  RS_SPC <- data.frame(RS = as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred),name = "Supervised principal components")
  rs <- cbind(mixed[,2:3],RS=as.numeric(superpc.predict(fit,data,list(x=t(mixed[,-c(1,2,3)]),y=mixed$OS.time,censoring.status=mixed$OS,featurenames=colnames(mixed)[-c(1,2,3)]),threshold = cv.fit$thresholds[which.max(cv.fit[["scor"]][1,])],n.components = 1)$v.pred))
  cc <- data.frame(name=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('Supervised PCA')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  
  
  #### plsRcox ####
  set.seed(seed)
  
  pdf('plsRcox.pdf')
  cv.plsRcox.res=cv.plsRcox(list(x=mixed[,-c(1,2,3)],time=mixed$OS.time,status=mixed$OS),nt=10,verbose = FALSE)
  fit <- plsRcox(mixed[,-c(1,2,3)],time=mixed$OS.time,event=mixed$OS,nt=as.numeric(cv.plsRcox.res[5]))
  RS_PLS <- data.frame(RS=as.numeric(predict(fit,type="lp",newdata=mixed[,-c(1,2,3)])),name='plsRcox')
  rs <- cbind(mixed[,2:3],RS=as.numeric(predict(fit,type="lp",newdata=mixed[,-c(1,2,3)])))
  cc <- data.frame(namx=as.numeric(summary(coxph(Surv(OS.time,OS)~RS,rs))$concordance[1]))
  colnames(cc) <- c('plsRcox')
  cc <- t(cc)
  colnames(cc) <- c('C')
  result <- rbind(result,cc)
  rm(list=c('fit'))
  dev.off()
  
}

