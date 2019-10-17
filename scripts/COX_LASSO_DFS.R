library(glmnet)
library(survival)
library(mise)
mise()
library(glmnet)
library(survival)
brca<-read.csv("./file.csv",header=T, na.strings=c("","NA"))
brca3<-brca[brca$DFS_MONTHS>0,]
brca4=brca3[!is.na(brca3$Sample_ID),]
brca4=brca4[!is.na(brca4$CLINICAL_STAGE),]
brca4$DFS_STATUS=as.character(brca4$DFS_STATUS)
brca4$DFS_STATUS[brca4$DFS_STATUS=="Recurred/Progressed"]<-"1"
brca4$DFS_STATUS[brca4$DFS_STATUS=="DiseaseFree"]<-"0"
brca4$DFS_STATUS=as.numeric(brca4$DFS_STATUS)
brca4<-na.omit(brca4)
y <- Surv(brca4$DFS_MONTHS,brca4$DFS_STATUS)
brca5<-brca4[,-c(83:86)]
brca5<-brca5[,-c(1:4)]
brca5$CLINICAL_STAGE=sub("*(A|B|C|a|b|c)","",brca5$CLINICAL_STAGE)
brca5$CLINICAL_STAGE=sub("*(A|B|C|a|b|c)","",brca5$CLINICAL_STAGE)
### scale the data frame
#brca6<-as.data.frame(scale(brca5))
brca6=brca5
#cox_formula <- as.formula(paste("~ ", paste(colnames(brca6), collapse="+")))
cox_formula <- as.formula(paste("~ ", paste(colnames(brca6), collapse="+"),"+0"))
cox_formula
brca6<-na.omit(brca6)
x <- model.matrix(cox_formula,brca6)
x[,c("AGE")]<-scale(x[,c("AGE")])
#cv.fit=cv.glmnet(x,y,family="cox", maxit = 1000)
fit=cv.glmnet(x,y,family="cox", nfolds=10, maxit = 1000)
important_features <- as.matrix(coef(fit, s = "lambda.min"))
##extract important features
brca7<-as.data.frame(x[,c("gain2q","gain10p","gain17q","gain21q","CLINICAL_STAGE")])
#brca7<-scale(brca7) ##scale only numeric part
colnames(brca7)<-gsub("CLINICAL_STAGEStge ","CLINICAL_STAGEStage_",colnames(brca7))
cox_formula <- as.formula(paste("y ~", paste(colnames(brca7), collapse="+")))
res.cox <- coxph(cox_formula, data =  brca7)
summary(res.cox)
##adjust for multiple testing
p.adjust(pval, method = "BH", n = length(res.cox$coefficients)) 
