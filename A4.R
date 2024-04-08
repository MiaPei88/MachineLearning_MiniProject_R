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
library(pROC)



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

##### Model Training #####

### Logistic Regression ###

# Creating train and test sets for logistic glmnet function
fit_train <- lm(population ~ ., data = train_set)
X_train <- model.matrix(fit_train)[,-1]
y_train <- train_set$population
fit_test <- lm(population ~ ., data = test_set)
X_test <- model.matrix(fit_test)[,-1]
y_test <- test_set$population

# Training logistic function via LASSO and measuring CV error with AUC
logistic_train <- cv.glmnet(X_train ,y_train, nfolds = 10, family="binomial", alpha=1, type.measure = "auc")
plot(logistic_train)

# Predicting on the test set using lambda.min and lambda.1se models
prds.test_min <- predict(logistic_train,newx = X_test, type = "response", s=logistic_train$lambda.min)[,1]
prds.test_1se <- predict(logistic_train,newx = X_test, type = "response", s=logistic_train$lambda.1se)[,1]

# Plotting ROC curves for lambda.min and lambda.1se test performance
par(mfrow=c(1,2))
auc.test_min <- roc(y_test,prds.test_min)
auc.test_min
plot(auc.test_min)
auc.test_1se <- roc(y_test,prds.test_1se)
auc.test_1se
plot(auc.test_min)

# Listing coefficients retained by lambda.min
coef.min <- coef(logistic_train,s=logistic_train$lambda.min)[,1]
coef.min[coef.min!=0]

### Random Forests ###
RF_train <- data.frame(lapply(train_set, as.factor))

RF_model <- randomForest(population ~ ., data = RF_train)
RF_model
plot(RF_model)

