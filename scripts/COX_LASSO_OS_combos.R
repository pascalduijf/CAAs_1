setwd("~/cal-scna-combos/OS_combos/")
library(mise)
mise()
library(glmnet)
library(survival)
brca<-read.csv("~/cal-scna-combos/",header=T, na.strings=c("","NA"))
brca3<-brca[brca$OS_MONTHS>0,]
brca4=brca3[!is.na(brca3$Sample_ID),]
brca4=brca4[!is.na(brca4$CLINICAL_STAGE),]
brca4$OS_STATUS=as.character(brca4$OS_STATUS)
brca4$OS_STATUS[brca4$OS_STATUS=="DECEASED"]<-"1"
brca4$OS_STATUS[brca4$OS_STATUS=="LIVING"]<-"0"
brca4$OS_STATUS=as.numeric(brca4$OS_STATUS)
brca4<-na.omit(brca4)
y <- Surv(brca4$OS_MONTHS,brca4$OS_STATUS)
brca5<-brca4[,-c(83:86)]
brca5<-brca5[,-c(1:4)]
brca5<-brca5[,-c(1:78)]
brca5$CLINICAL_STAGE=sub("*(A|B|C|a|b|c)","",brca5$CLINICAL_STAGE)
brca5$CLINICAL_STAGE=sub("*(A|B|C|a|b|c)","",brca5$CLINICAL_STAGE)
brca6<-brca5
#cox_formula <- as.formula(paste("~ ", paste(colnames(brca6), collapse="+")))
cox_formula <- as.formula(paste("~ ", paste(colnames(brca6), collapse="+"),"+0"))
cox_formula
brca6<-na.omit(brca6)
x <- model.matrix(cox_formula,brca6)
x[,c("AGE")]<-scale(x[,c("AGE")])
cv.fit=cv.glmnet(x,y,family="cox", maxit = 10000)
#fit=SIS(x, y, family = c("cox"),penalty = c("SCAD", "MCP", "lasso"), concavity.parameter = switch(penalty,SCAD = 3.7, 3), tune = c("bic", "ebic", "aic", "cv"), nfolds = 10,type.measure = c("deviance", "class", "auc", "mse", "mae"),gamma.ebic = 1, nsis = NULL, iter = TRUE,iter.max = ifelse(greedy ==FALSE, 10, rfloor(nrow(x)/log(nrow(x)))), varISIS = c("vanilla", "aggr","cons"), perm = FALSE, q = 1, greedy = FALSE,greedy.size = 1,seed = 0, standardize = TRUE)

#fit=cv.glmnet(x,y,family="cox", maxit = 100000000)
important_features <- as.matrix(coef(cv.fit, s = "lambda.min"))
View(important_features)
important_features1=as.data.frame(important_features)
important_features1$combos=row.names(important_features1)
colnames(important_features1)=c("importance","combos")
#row.names(important_features1)=NULL
important_features_not_zero=as.data.frame(important_features1$importance >0 | important_features1$importance <0,important_features1$combos)
View(important_features_not_zero)
important_features_not_zero$importance_category=as.character(important_features_not_zero$importance_category)
important_features_not_zero$combos=rownames(important_features_not_zero)
rownames(important_features_not_zero)=NULL
important_features_not_zero_=as.data.frame(important_features_not_zero[important_features_not_zero$importance_category=="TRUE",])
important_features_not_zero_
colnames(important_features_not_zero)=c("importance","combos")
important_features_not_zero$importance=as.character(important_features_not_zero$importance)
important_features_not_zero_=as.data.frame(important_features_not_zero[important_features_not_zero$importance=="TRUE",])
##extract important features
brca7<-as.data.frame(x[,important_features_not_zero_$combos])
#brca7<-scale(brca7) ##scale only numeric part
colnames(brca7)<-gsub("CLINICAL_STAGEStge ","CLINICAL_STAGEStage_",colnames(brca7))
###if only one column is present then
#colnames(brca7)<-"the name of that one covariate"
cox_formula <- as.formula(paste("y ~", paste(colnames(brca7), collapse="+")))
res.cox <- coxph(cox_formula, data =  brca7)
summary(res.cox)
##adjust for multiple testing
p.adjust(pval, method = "BH", n = length(res.cox$coefficients)) 
pu=publish(res.cox,org=T)
poka=as.data.frame(pu$regressionTable)
poka2<-subset(poka,poka$`p-value`<.05)
publish(poka2,org=T)

