# R script for determining synthetic lethal, 
# synergistic resistance and co-lost and 
# co-gained focal copy number alterations on 
# the same chromosome arm
table <- read.table("~/main_table.csv", header=TRUE, sep=",", check.names = FALSE)
table <- table[-c(1:2,4:5)]
table1 <- subset(table, TCGA_CODE == "ALL")
table1 <- table1[-1]
table2 <- subset(table, TCGA_CODE == "BLCA")
table2 <- table2[-1]
table3 <- subset(table, TCGA_CODE == "BRCA")
table3 <- table3[-1]
table4 <- subset(table, TCGA_CODE == "CESC")
table4 <- table4[-1]
table5 <- subset(table, TCGA_CODE == "COADREAD")
table5 <- table5[-1]
table6 <- subset(table, TCGA_CODE == "DLBC")
table6 <- table6[-1]
table7 <- subset(table, TCGA_CODE == "ESCA")
table7 <- table7[-1]
table8 <- subset(table, TCGA_CODE == "GBM")
table8 <- table8[-1]
table9 <- subset(table, TCGA_CODE == "HNSC")
table9 <- table9[-1]
table10 <- subset(table, TCGA_CODE == "KIRC")
table10 <- table10[-1]
table11 <- subset(table, TCGA_CODE == "LAML")
table11 <- table11[-1]
table12 <- subset(table, TCGA_CODE == "LCML")
table12 <- table12[-1]
table13 <- subset(table, TCGA_CODE == "LGG")
table13 <- table13[-1]
table14 <- subset(table, TCGA_CODE == "LIHC")
table14 <- table14[-1]
table15 <- subset(table, TCGA_CODE == "LUAD")
table15 <- table15[-1]
table16 <- subset(table, TCGA_CODE == "LUSC")
table16 <- table16[-1]
table17 <- subset(table, TCGA_CODE == "MB")
table17 <- table17[-1]
table18 <- subset(table, TCGA_CODE == "MESO")
table18 <- table18[-1]
table19 <- subset(table, TCGA_CODE == "MM")
table19 <- table19[-1]
table20 <- subset(table, TCGA_CODE == "NB")
table20 <- table20[-1]
table21 <- subset(table, TCGA_CODE == "OV")
table21 <- table21[-1]
table22 <- subset(table, TCGA_CODE == "PAAD")
table22 <- table22[-1]
table23 <- subset(table, TCGA_CODE == "PRAD")
table23 <- table23[-1]
table24 <- subset(table, TCGA_CODE == "SCLC")
table24 <- table24[-1]
table25 <- subset(table, TCGA_CODE == "SKCM")
table25 <- table25[-1]
table26 <- subset(table, TCGA_CODE == "STAD")
table26 <- table26[-1]
table27 <- subset(table, TCGA_CODE == "THCA")
table27 <- table27[-1]
table28 <- subset(table, TCGA_CODE == "UCEC")
table28 <- table28[-1]
table29 <- table[-1]
table.function=function(t1,t2) {
    matrix1=matrix(NA,2,2)
    matrix1[1,1]=sum(t1==0 & t2==0)
    matrix1[1,2]=sum(t1==0 & t2==1)
    matrix1[2,1]=sum(t1==1 & t2==0)
    matrix1[2,2]=sum(t1==1 & t2==1)
    return(matrix1)
}
# Code below was run for table1-table28; sample code below only shows code for table4 as an example
matrix1=matrix(NA,265*choose(788,2),8) # no of drugs, no of features, choose combos of 2
colnames(matrix1)=c("feature_1","feature_2","drug","R_both_features",
                "R_not_both_features","S_both_features",
                "S_not_both_features","fisher.test.pval")
m=1
for (k in 1:265) { #drugs
  for (i in 1:755) { #first feature of pair
    for (j in (i+1):755) { #second feature of pair
      tmp1=which(table4[,k+755]=="R") #vector with row numbers that are "R"
      tmp2=which(table4[,k+755]=="S") #vector with row numbers that are "S"
      tmp3=table.function(table4[tmp1,i],table4[tmp1,j]) #table for "R"
      tmp4=table.function(table4[tmp2,i],table4[tmp2,j]) #table for "S"
      tmp5=matrix(c(tmp3[2,2],sum(tmp3)-tmp3[2,2],
                    tmp4[2,2],sum(tmp4)-tmp4[2,2]),nrow=2,byrow=T) #2x2 matrix
      matrix1[m,]=c(colnames(table4)[i],
                colnames(table4)[j],
                colnames(table4)[k+755],
                tmp5[1,1],tmp5[1,2],tmp5[2,1],tmp5[2,2],
                fisher.test(tmp5)$p.value)
      m=m+1
    }
  }
}
matrix1 <- na.omit(matrix1[matrix1[, "fisher.test.pval"] < 0.05,])
write.csv(matrix1, file = "~/Synth_Lethal_table4.csv", row.names=F)
# Files "Synth_Lethal_table1.csv" - "Synth_Lethal_table28.csv" were combined for Supplemental Table