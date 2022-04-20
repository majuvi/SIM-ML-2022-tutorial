library(glmnet)
library(dplyr)
library(DALEX)
library(pROC)
library(ranger)

set.seed(321)

# load data
data(titanic_imputed, package = "DALEX")
data_titanic <- titanic_imputed

# analyse data
dim(data_titanic)
head(data_titanic)
summary(data_titanic)
str(data_titanic)


# scale data
columns_scale = c("age", "fare", "sibsp", "parch")
data_titanic[, columns_scale] <- scale(data_titanic[, columns_scale])


# split data in train and validation set
index_train_test_split = sample(nrow(data_titanic), 0.8*nrow(data_titanic))

titanic_train <- data_titanic[index_train_test_split,]
titanic_val <- data_titanic[-index_train_test_split,]

X_train <- dplyr::select(titanic_train, -survived)
y_train <- titanic_train$survived

X_val <- dplyr::select(titanic_val, -survived)
y_val <- titanic_val$survived


# train model
fit <- glmnet(X_train, y_train, family="binomial", lambda=0.01, alpha=0.5)

# evaluation fit 
coef(fit)

y_true = titanic_val$survived  
y_pred = as.factor(predict(fit, data.matrix(X_val), type="class"))
        
cm <- table(y_pred, y_true)
print(cm)

accuracy <- sum(cm[1], cm[4]) / sum(cm[1:4])
print(accuracy)

# calculating ROC and AUC
y_pred_prob = as.vector(predict(fit, data.matrix(X_val), type="response"))

proc <- pROC::roc(y_true, y_pred_prob)
plot(proc)
pROC::auc(y_true, y_pred_prob)


# take whole dataset for cross validation
titanic_X <- data.matrix(dplyr::select(data_titanic, -survived))
titanic_y <- titanic_imputed$survived

# cross validation
for(alpha_cv in c(0, 0.5, 1)){
  fit_cv <- cv.glmnet(titanic_X, titanic_y, nfold=5, family="binomial", type.measure = "auc", alpha=alpha_cv)
  
  plot(fit_cv, main=paste("alpha =", toString(alpha_cv)))

  max_auc <- fit_cv$cvm[fit_cv$lambda == fit_cv$lambda.min]
  
  print(paste("alpha=", alpha_cv, "lambda_min=", fit_cv$lambda.min, "AUC=", max_auc)) # print MSE
  
}



#### random forest ####
# reload data
data(titanic_imputed, package = "DALEX")
data_titanic <- titanic_imputed

# split to train and validation set 
index_train_test_split = sample(nrow(data_titanic), 0.8*nrow(data_titanic))
titanic_train <- data_titanic[index_train_test_split,]
titanic_val <- data_titanic[-index_train_test_split,]

# fit random forest
fit_rf <- ranger::ranger(survived ~ ., data = titanic_train, 
                         classification = TRUE, num.trees = 50, 
                         probability = TRUE)

# make predictions
pred_rf <- predict(fit_rf, data = titanic_val)
y_pred <- pred_rf$predictions
y_true <- titanic_val$survived

# ROC and AUC
proc <- pROC::roc(y_true, y_pred[,2])
plot(proc)
pROC::auc(y_true, y_pred[,2])


# evaluate
cm <- table(y_true, y_pred[,2] > 0.5)
accuracy <- sum(cm[1], cm[4]) / sum(cm[1:4])
print(accuracy)








