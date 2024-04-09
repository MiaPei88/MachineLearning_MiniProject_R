##### Data Preparation #####
## Install and load required packages
#install.packages("vcfR")
library(vcfR)
library(tidyverse)
library(dplyr)
library(stringr)
library(randomForest)
library(ggplot2)
library(glmnet)
library(caret)
library(pROC)
library(reshape2)
library(cvms)
library(rpart)
library(rpart.plot)
library(rattle)

## Read in all the vcf files and extract the genotype information
vcf_LCT_EUR <- read.vcfR("./filtered_LCT_EUR.vcf")
geno_LCT_EUR <- extract.gt(vcf_LCT_EUR, element = "GT")

vcf_LCT_EAS <- read.vcfR("./filtered_LCT_EAS.vcf")
geno_LCT_EAS <- extract.gt(vcf_LCT_EAS, element = "GT")

vcf_MCM6_EUR <- read.vcfR("./filtered_MCM6_EUR.vcf")
geno_MCM6_EUR <- extract.gt(vcf_MCM6_EUR, element = "GT")

vcf_MCM6_EAS <- read.vcfR("./filtered_MCM6_EAS.vcf")
geno_MCM6_EAS <- extract.gt(vcf_MCM6_EAS, element = "GT")

## Write a function to convert the matrix into a data frame with proper format
convert_geno <- function(geno, pop, gene) {
  geno <- as.data.frame(geno, stringsAsFactors = FALSE)
  
  # Change all the values in the data frame into proper format
  geno <- geno %>%
    mutate(across(everything(), ~case_when(
    . == "0|0" ~ 0,
    . == "0|1" | . == "1|0" ~ 1,
    . == "1|1" ~ 2,
    TRUE ~ NA_real_)))
  
  # Transpose individuals as rows and SNPs as columns
  geno_transposed <- as.data.frame(t(geno))
  
  # Add gene names as prefix of column names
  colnames(geno_transposed) <- paste(gene, colnames(geno_transposed), sep = "_")
  
  # Add two columns, population and individual
  geno_transposed <- geno_transposed %>%
    mutate(population = pop) %>%
    mutate(individual = row.names(.))
  
  return(geno_transposed)
}

## Convert the 4 matrices into data frames with appropriate format
converted_LCT_EUR <- convert_geno(geno_LCT_EUR, "EUR", "LCT")
converted_LCT_EAS <- convert_geno(geno_LCT_EAS, "EAS", "LCT")
converted_MCM6_EUR <- convert_geno(geno_MCM6_EUR, "EUR", "MCM6")
converted_MCM6_EAS <- convert_geno(geno_MCM6_EAS, "EAS", "MCM6")

## Use full_join to combine data frames from the same population by their common columns
common_cols_EUR <- intersect(colnames(converted_LCT_EUR), colnames(converted_MCM6_EUR))
common_cols_EAS <- intersect(colnames(converted_LCT_EAS), colnames(converted_MCM6_EAS))
combined_EUR <- full_join(converted_LCT_EUR, converted_MCM6_EUR, by = common_cols_EUR)
combined_EAS <- full_join(converted_LCT_EAS, converted_MCM6_EAS, by = common_cols_EAS)

## Combine the two data frames of EUR and EAS together by their common columns
common_cols <- intersect(colnames(combined_EUR), colnames(combined_EAS))
combined_df <- full_join(combined_EUR, combined_EAS, by = common_cols)
# Change the row names to individual id
rownames(combined_df) <- combined_df$individual

## Remove columns with NA values
combined_df_commonSNP <- combined_df %>%
  # remove "individual" column
  select(-individual) %>%
  select(where(~ !any(is.na(.)))) %>%
  # move "population" column to the front
  select(population, everything())

combined_logit <- combined_df %>%
  # remove "individual" column
  select(-individual) %>%
  select(where(~ !any(is.na(.)))) %>%
  # move "population" column to the front
  select(population, everything()) %>%
  mutate(population = case_when(
    population == "EUR" ~ 0,
    population == "EAS" ~ 1,
    TRUE ~ as.integer(NA)
  ))


## Data splitting
set.seed(3575) 
train.index <- sample(1:nrow(combined_df_commonSNP), round(0.70*nrow(combined_df_commonSNP)))
train_set <- combined_df_commonSNP[train.index,]
test_set <- combined_df_commonSNP[-train.index,]
train_set_logit <- combined_logit[train.index,]
test_set_logit <- combined_logit[-train.index,]

##### Model Training #####

### Logistic Regression ###

# Creating train and test sets for logistic glmnet function
fit_train <- lm(population ~ ., data = train_set_logit)
X_train <- model.matrix(fit_train)[,-1]
y_train <- train_set_logit$population
fit_test <- lm(population ~ ., data = test_set_logit)
X_test <- model.matrix(fit_test)[,-1]
y_test <- test_set_logit$population

# Training logistic function via LASSO and measuring CV error with AUC
set.seed(3575)
logistic_train <- cv.glmnet(X_train ,y_train,
                            nfolds = 10, family="binomial",
                            alpha=1, type.measure = "auc")
par(mfrow=c(1,1))
plot(logistic_train)

# Predicting on the test set using lambda.min and lambda.1se models
prds.test_min <- predict(logistic_train,
                         newx = X_test,
                         type = "response",
                         s=logistic_train$lambda.min)[,1]
prds.test_1se <- predict(logistic_train,
                         newx = X_test,
                         type = "response", 
                         s=logistic_train$lambda.1se)[,1]

# Plotting ROC curves for lambda.min and lambda.1se test performance
par(mfrow=c(1,2))
auc.test_min <- roc(y_test,prds.test_min)
auc.test_min
plot(auc.test_min)
snsp.test_min <- cbind(auc.test_min$sensitivities,auc.test_min$specificities)
indx <- which.max(apply(snsp.test_min,1,min))
abline(h=snsp.test_min[indx,1],v=snsp.test_min[indx,2], col='blue', lty=2)

auc.test_1se <- roc(y_test,prds.test_1se)
auc.test_1se
plot(auc.test_1se)
snsp.test_1se <- cbind(auc.test_1se$sensitivities,auc.test_1se$specificities)
indx2 <- which.max(apply(snsp.test_1se,1,min))
abline(h=snsp.test_1se[indx2,1],v=snsp.test_1se[indx2,2], col='blue', lty=2)
par(mfrow=c(1,1))


# Listing coefficients retained by lambda.min
coef.min <- coef(logistic_train,s=logistic_train$lambda.min)[,1]
coef.min[coef.min!=0]
sort(abs(coef.min[coef.min!=0]), decreasing = T)
names(sort(abs(coef.min[coef.min!=0]), decreasing = T))[-1]

#Confusion Matrix for lambda.min
pred_class_lasso <- ifelse(prds.test_min > 0.5, 1, 0)
matrix_table_lasso <- table(observed = y_test,
                      predicted = pred_class_lasso)

# Store the table as a tibble
matrix_tibble_lasso <- as_tibble(matrix_table_lasso)

# Plot the tibble as a confusion matrix plot
cm_lasso <- plot_confusion_matrix(matrix_tibble_lasso,
                                  target_col = "observed",
                                  prediction_col = "predicted",
                                  counts_col = "n",
                                  add_sums = TRUE,
                                  add_normalized = FALSE,
                                  add_col_percentages = FALSE,
                                  add_row_percentages = FALSE,
                                  palette = "Purples",
                                  sums_settings = sum_tile_settings(palette = "Oranges",
                                                                    label = "Total"))

plot(cm_lasso)

### Random Forests ###
## Train the model
train_set$population <- as.factor(train_set$population)
test_set$population <- as.factor(test_set$population)
set.seed(3575) 
RF_model <- randomForest(population ~ ., data = train_set, keep.forest = TRUE)
forest_id <- RF_model$forest

# Check the model training result
RF_model

# Draw the OOB error rate plot of the model
plot(RF_model)
legend("topright", 
       # Use the class names from the model to specify the colors of line
       legend=colnames(RF_model$err.rate),  
       col=c("black","red", "green"), 
       lty=1, 
       cex=0.6)

# Check the variable importance
varImpPlot(RF_model)


## Tune parameter mtry 
set.seed(3575)
trcontrol <- trainControl(method='repeatedcv', 
                          number=10, 
                          repeats=1, 
                          summaryFunction = twoClassSummary, 
                          classProbs = TRUE,
                          search='grid')

# Check the square root of number of variables and set the tuning range of mtry
round(sqrt(ncol(train_set)-1))
tunegrid <- expand.grid(mtry = c(1:11)) 

# Train the model using different mtry until it finds the best
system.time(tuned_RF <- train(population ~ ., 
                  data = train_set,
                  trControl = trcontrol,
                  method = "rf",
                  metric = "ROC",
                  tuneGrid = tunegrid))

# Check the tuned model
print(tuned_RF)					   
plot(tuned_RF)


## Prediction
# Predict on test set
pred_test_rf <- predict(tuned_RF, newdata = test_set, type="prob")[,1]

# Compute the confusion matrix
pred_class_rf <- ifelse(pred_test_rf > 0.5, "EAS", "EUR")
conf_mat <- confusionMatrix(as.factor(pred_class_rf), test_set$population)
conf_mat

# Visualize the confusion matrix
conf_mat_melted <- as_tibble(conf_mat$table)
plot_confusion_matrix(conf_mat_melted, target_col = "Reference", 
                      prediction_col = "Prediction", counts_col = "n", 
                      add_sums = TRUE, add_normalized = FALSE,
                      add_col_percentages = FALSE,
                      add_row_percentages = FALSE, palette = "Purples", 
                      sums_settings = sum_tile_settings(palette = "Oranges", 
                                                        label = "Total"))


## ROC-AUC
# Calculate ROC of the prediction on test set
roc_rf <- roc(test_set$population, pred_test_rf)

# Sensitivity and Specificities
snsp_rf <- cbind(roc_rf$sensitivities, roc_rf$specificities)
indx_rf <- which.max(apply(snsp_rf,1,min))
cutoff_rf <- roc_rf$thresholds[indx_rf]
indx_rf
snsp_rf[indx_rf,]

# Visualize the ROC-AUC plot
plot(roc_rf, main = "ROC Curve")
abline(h=snsp_rf[indx_rf,1],v=snsp_rf[indx_rf,2], col='blue', lty=2)
text(0.6, 0.6, paste("AUC =", round(auc(roc_rf), 3)), col = "#ff0000")


### Decision Tree ###
## Model Training ##
dt_model <- rpart(population ~ ., data = train_set, method = "class")
par(mfrow = c(1,2), cex = 1)
summary(dt_model)
fancyRpartPlot(dt_model, main = "Decision Tree Visualization for Full Tree")


## Tree Pruning ##
prune.cp <- function(cptable){
  CP <- cptable[,1]
  cp <- sqrt(CP * c(Inf, CP[-length(CP)]))
  xerr <- cptable[,4]
  minrow <- which.min(xerr)
  
  xerr.1se <- sum(cptable[minrow,4:5])
  index <- min(which(xerr < xerr.1se))
  
  cp.1se <- cp[index]
  return(as.numeric(cp.1se) )}

prune.cp(dt_model$cptable)
pruned_tree <- prune(dt_model, cp = prune.cp(dt_model$cptable))
summary(pruned_tree)
par(mfrow = c(1,1), cex = 1.2)
fancyRpartPlot(pruned_tree, main = "Decision Tree Visualization for Pruned Tree")


## Model Validation ##
# Predict using the pruned decision tree model on the test set
pred_dt <- predict(pruned_tree, newdata = test_set, type = "class")


## Confusion Matrix Plot
# Create a confusion matrix in the table format
matrix_table <- table(observed = test_set$population, predicted = pred_dt)

# Store the table as a tibble
matrix_tibble <- as_tibble(matrix_table)

# Plot the tibble as a confusion matrix plot
par(mfrow = c(1,2), cex = 1)
cm_plot_3 <- plot_confusion_matrix(matrix_tibble, target_col = "observed",
                                   prediction_col = "predicted", counts_col = "n",
                                   add_sums = TRUE, add_normalized = FALSE,
                                   add_col_percentages = FALSE, add_counts = TRUE,
                                   add_row_percentages = FALSE, palette = "Purples",
                                   sums_settings = sum_tile_settings(palette = "Oranges",
                                                                     label = "Total")) +
  ggtitle("Confusion Matrix Plot for Decision Tree Classifier")
cm_plot_3


## ROC-AUC
# Calculate ROC of the prediction on test set
roc_dt <- roc(test_set$population, as.numeric(pred_dt))

# Sensitivity and Specificities
snsp_dt <- cbind(roc_dt$sensitivities, roc_dt$specificities)
indx_dt <- which.max(apply(snsp_dt,1,min))
cutoff_dt <- roc_dt$thresholds[indx_dt]
indx_dt
snsp_dt[indx_dt,]

# Visualize the ROC-AUC plot
par(mfrow = c(1,1), cex = 1)
plot(roc_dt, main = "ROC Curve for Decision Tree Classifier")
abline(h = snsp_dt[indx_dt,1], v = snsp_dt[indx_dt,2], col = 'blue', lty = 2)
text(0.6, 0.6, paste("AUC =", round(auc(roc_dt), 3)), col = "#ff0000")

