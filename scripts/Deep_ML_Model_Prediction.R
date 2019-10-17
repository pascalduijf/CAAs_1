suppressMessages(library(yardstick)) # used to check model performance
suppressMessages(library(tidyverse)) # implement ggplot functionality for beautiful and meaningful visualization
suppressMessages(library(caret))
setwd("./")
suppressPackageStartupMessages({ library (readxl)
  library (keras);library (lime);library (dplyr)
  ;library (recipes);library (yardstick)
  ;library(knitr);library (DT) })
#remove saved model generated from previous run  
system("rm best_model.h5")
ca<-read.csv("data_CALSCNAsONLY_          imputed_39.csv")
ca=ca[,-1]
names(ca)
cal_erlotinib<-ca[,1:90]
bad=sapply(cal_erlotinib, function(x) all(is.na(x)))
cal_erlotinib<-cal_erlotinib[,!bad]
names(cal_erlotinib)		  
cal_erlotinib<-cal_erlotinib[,-89]
names(cal_erlotinib)		  
# Spliting the data
set.seed(100) # For Reproducibility
train <- sample(1:nrow(cal_erlotinib),nrow(cal_erlotinib)*.8) # sample of the specified size from the total row
#train<-train[,!bad] ##remove na features
model_data <- cal_erlotinib[train,]
test_data <- cal_erlotinib[-train,]
# Predictor variable for the training and testing set
X_train <- model_data[,c(1:88)]
X_test <- test_data[,c(1:88)]
##Scaling the predictor valriable
#X_train <-scale(X_train)
# Target Variable for training and testing set
y_train <- model_data[,-c(1:88)]
y_test <- test_data[,-c(1:88)]
build_model <- function() {
    model <- keras_model_sequential() %>%
        layer_dense(units = 256,activation = "relu",
                    input_shape = dim(X_train)[2]) %>%  layer_dropout(rate = 0.5) %>%
        layer_dense(units = 132, activation = "relu") %>%  layer_dropout(rate = 0.4) %>%
        layer_dense(units = 64, activation = "relu") %>%  layer_dropout(rate = 0.3) %>%
        layer_dense(units = 32, activation = "relu") %>%  layer_dropout(rate = 0.2) %>%  
        layer_dense(units = 16, activation = "relu") %>%  layer_dropout(rate = 0.1) %>%
        layer_dense(units = 8, activation = "relu") %>%  layer_dense(units = 1,activation="sigmoid")
    model %>% compile(
        loss = "binary_crossentropy",
        optimizer = optimizer_rmsprop(lr=.0001),
        metrics = list("accuracy")
    )
    model
}
model <- build_model()
model %>% summary()
mc = callback_model_checkpoint('best_model.h5', monitor='val_acc', mode='max', verbose=1, save_best_only = T)
k=callback_early_stopping(monitor='val_acc',verbose=1,mode = "max",patience=300)
history <- fit (
    object           = model,             # => Our Model
    x                = as.matrix (X_train), #=> Matrix
    y                = y_train,             #=> Numeric Vector
    batch_size       = 555,     #=> #OfSamples/gradient update in each epoch
    epochs           = 300,     #=> Control Training cycles
    validation_split = 0.2, callbacks = c(mc,k))  #=>
model = load_model_hdf5('best_model_caa.h5')
yhat_keras_class_vec <-
predict_classes (object = model,
x = as.matrix(X_test)) %>%
as.vector()
yhat_keras_prob_vec <-
predict_proba (object = model,
x = as.matrix(X_test)) %>%
as.vector()
estimates_keras_tbl <- tibble(
truth      = as.factor(y_test) %>%
fct_recode (Resistance = "1", Sensitive = "0"),
estimate   = as.factor(yhat_keras_class_vec) %>%
fct_recode (Resistance = "1", Sensitive = "0"),
class_prob = yhat_keras_prob_vec)
options(scipen = 999)
head(estimates_keras_tbl, 20)
confusion<-confusionMatrix(estimates_keras_tbl$estimate, estimates_keras_tbl$truth, positive = "Sensitive")
recall <- sensitivity(estimates_keras_tbl$estimate, estimates_keras_tbl$truth, positive="Sensitive")
precision <-posPredValue(estimates_keras_tbl$estimate, estimates_keras_tbl$truth, positive="Sensitive")
F1 <- (2 * precision * recall) / (precision + recall)
file <- file("./confusion_matrices/confusion_matrix_F1_Recall_Precision_caa_drug.txt")				 
writeLines(F1, file)
close(file)