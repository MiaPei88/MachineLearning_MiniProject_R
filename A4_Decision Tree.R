##### Data Preparation #####
## Install and load required packages
library(rpart)
library(rpart.plot)
library(rattle)



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
pruned_tree <- prune(dt_model, cp = prune.cp(tree_model$cptable))
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

