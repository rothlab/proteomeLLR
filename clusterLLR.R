library(mclust)
library(ggbiplot)
library(umap)
library(naivebayes)

# This script takes a csv with LLRs at 0.1 VARITY score intervals for a set of genes 
# It will attempt to cluster the genes based on their LLR curves and output PCA 
# and UMAP plots
# Additionally, it attempts to build a classification model using naive bayes to 
# see if VARITY performs better for certain subgroups of genes (ex. inheritance, 
# different functional sub classes of genes)
#
# See sample input files under demo_data

# Read in LLRs at 0.1 score intervals for every gene with sufficient reference variants
llrs <- read.csv("LLR_intervals.csv", header = TRUE, check.names = FALSE)

# Look at the PCA plot without group classifications
llrs.pca <- prcomp(llrs[,2:ncol(llrs)], scale=F)
print(ggbiplot(llrs.pca, var.axes = F, choices = 1:2))

# Perform clustering with n groups (where n ranges from 1 to 6 in this case)
# and then visualize PCA plots and UMAPs with n groups
my_colours <- c("red", "blue", "green", "black", "purple", "brown")

for (i in 1:6) {
  llrs.clust <- Mclust(llrs[,2:ncol(llrs)], G=i)
  
  llrs["classification"] = llrs.clust$classification
  classification_colors = my_colours[llrs$classification]
  
  # Visualize PCA plot with group classifications
  p <- ggbiplot(llrs.pca, var.axes = F, choices = 1:2, groups = factor(llrs$classification))
  pdf(paste("llrs_pca_", i, "_groups.pdf", sep = ""),15,5)
  op <- par(mfrow=c(1,1))
  print(p)
  par(op)
  dev.off()
  
  # Visualize UMAP with group classifications
  pdf(paste("llrs_umap_", i, "_groups.pdf", sep = ""),15,5)
  op <- par(mfrow=c(1,1))
  llrs.umap <- umap(llrs[,2:ncol(llrs)])
  plot(llrs.umap$layout,col=my_colours[llrs$classification],pch=20, 
       main = "UMAP visualization of LLR dataset", xlab="", ylab="")
  par(op)
  dev.off()
}

#' Return the gene closest to the centroid of a cluster
#' 
#' @param llrs dataframe with llr data
#' @param clust cluster object
#' @param group the group of interest within clust
#' 
#' @return a vector with the distance between the centroid and the closest gene 
#' and the gene name
#' @export
#' 
findClusterCenter <- function(llrs, clust, group) {
  minDiffGene = ""
  minDiff = Inf
  clusterMean = clust$parameters$mean[, group] 
  
  for (i in 1:nrow(llrs[llrs$classification == group,])) {
    geneName = llrs[i,1]
    geneLLRs = llrs[i,2:12]
    
    diff = sum(sqrt((geneLLRs - clusterMean)^2))
    
    if (diff < minDiff) {
      minDiffGene = geneName
      minDiff = diff
    }
  }
  return(c(minDiff, minDiffGene))
}

################ Classify genes based on inheritance ###########################

# Read in LLRs at 0.1 score intervals for every gene with sufficient reference 
# variants WITH inheritance annotations from GenCC
inheritanceLLR <- read.csv("LLR_intervals_w_inheritance.csv", header = TRUE, check.names = FALSE)
inheritanceLLR$Inheritance <- as.factor(inheritanceLLR$Inheritance)

ar = inheritanceLLR[inheritanceLLR$Inheritance == "AR",]
ad = inheritanceLLR[inheritanceLLR$Inheritance == "AD",]

# Get 'signature' LLR curve for each mode of inheritance
arLLRSignature = colMeans(ar[,2:12])
adLLRSignature = colMeans(ad[,2:12])
plot(arLLRSignature, adLLRSignature, main="'Signature' LLR curves of AR vs AD genes") # the signatures for AR and AD genes are similar

# Split data into train and test sets
split <- sample(2, nrow(inheritanceLLR), replace = T, prob = c(0.8, 0.2))
train <- inheritanceLLR[split == 1, 2:ncol(inheritanceLLR)]
test <- inheritanceLLR[split == 2, 2:ncol(inheritanceLLR)]

# Use 10-fold cross-validation from training set 
classifiedRate = c()
trainSplit <- sample(10, nrow(train), replace = T, prob = rep(c(0.1), 10))
for (i in 1:10) {
  trainSubset <- train[trainSplit != i, 2:ncol(train)]
  testSubset <- train[trainSplit == i, 2:ncol(train)]
  
  model <- naive_bayes(Inheritance ~ ., data = trainSubset, usekernel = T)
  predictions <- predict(model, testSubset[,-c(ncol(testSubset))])
  table1 <- table(predictions, testSubset$Inheritance)

  classifiedRate = c(classifiedRate, sum(diag(table1))/sum(table1))
}
sprintf("The average correct classification rate: %f", mean(classifiedRate))

# Compare performance with random predictions
randomClassifiedRate = c()
for (i in 1:50){
  randomPredictions = sample(c("AD", "AR", "other"), nrow(train), replace=T, 
                           prob=c(nrow(train[train$Inheritance == "AR",])/nrow(train), 
                                  nrow(train[train$Inheritance == "AD",])/nrow(train), 
                                  nrow(train[train$Inheritance == "other",])/nrow(train)))
  actual = as.vector(train$Inheritance)
  randomPredictionsTable = table(actual, randomPredictions)
  randomClassifiedRate = c(randomClassifiedRate, sum(diag(randomPredictionsTable))/sum(randomPredictionsTable))
}
sprintf("The average classification rate with random predictions on training data: %f", mean(randomClassifiedRate))

# The performance of our model is not that much better than random predictions...

# Try out model performance with test set
model <- naive_bayes(Inheritance ~ ., data = train, usekernel = T)
testPredictions <- predict(model, test[,-c(ncol(test))])
inheritanceConfusionMatrix <- table(testPredictions, test$Inheritance)
testClassifiedRate = sum(diag(inheritanceConfusionMatrix))/sum(inheritanceConfusionMatrix)
sprintf("The correct classification rate on the test data: %f", testClassifiedRate)


############### Classify genes based on functional activity ####################

# Read in LLRs at 0.1 score intervals for every gene with sufficient reference 
# variants WITH functional annotations from DAVID
annotatedLLR = read.csv("LLR_intervals_w_david_functional_annotation.csv", header = TRUE, check.names = FALSE)

# Check which functional activity groups have less than 10 genes - exclude these groups
groups = summary(as.factor(annotatedLLR$functional_cluster))
groups = groups[groups >= 10]
filteredAnnotatedLLR = annotatedLLR[annotatedLLR$functional_cluster %in% names(groups),]
filteredAnnotatedLLR$functional_cluster <- as.factor(filteredAnnotatedLLR$functional_cluster)

# Get 'signature' LLR curve for each functional activity cluster (with at least 10 genes)
cluster1 = filteredAnnotatedLLR[filteredAnnotatedLLR$functional_cluster == 1,] # ion channel activity
cluster2 = filteredAnnotatedLLR[filteredAnnotatedLLR$functional_cluster == 2,] # protein kinase activity
cluster6 = filteredAnnotatedLLR[filteredAnnotatedLLR$functional_cluster == 6,] # iron binding
cluster7 = filteredAnnotatedLLR[filteredAnnotatedLLR$functional_cluster == 7,] # neurotransmitter receptor activity
cluster12 = filteredAnnotatedLLR[filteredAnnotatedLLR$functional_cluster == 12,] # DNA-binding

cluster1_fingerprint = colMeans(cluster1[,2:12])
cluster2_fingerprint = colMeans(cluster2[,2:12])
cluster6_fingerprint = colMeans(cluster6[,2:12])
cluster7_fingerprint = colMeans(cluster7[,2:12])
cluster12_fingerprint = colMeans(cluster12[,2:12])
plot(x=names(cluster1_fingerprint), y=cluster1_fingerprint, xlab = "VARITY score bin", ylab = "LLR", type="l", main="Signature LLR curves of functionally-clustered genes") 
lines(x=names(cluster2_fingerprint), y=cluster2_fingerprint, col = "pink") 
lines(x=names(cluster6_fingerprint), y=cluster6_fingerprint, col = "red")
lines(x=names(cluster7_fingerprint), y=cluster7_fingerprint, col = "blue")
lines(x=names(cluster12_fingerprint), y=cluster12_fingerprint, col = "green")
legend(x="topleft",legend=c("Ion Channel Activity", "Protein Kinase Activity", "Iron Binding", "Neurotransmitter Receptor Activity", "DNA-Binding"),
       col=c("black", "pink", "red", "blue", "green"), lty=1, cex=0.8)
# the signatures for the 5 functional activity clusters are similar

# Split data into training and testing sets
split <- sample(2, nrow(filteredAnnotatedLLR), replace = T, prob = c(0.8, 0.2))
train <- filteredAnnotatedLLR[split == 1, 2:ncol(filteredAnnotatedLLR)]
test <- filteredAnnotatedLLR[split == 2, 2:ncol(filteredAnnotatedLLR)]

# Using 10-fold cross-validation with training set
classifiedRate = c()
trainSplit <- sample(10, nrow(train), replace = T, prob = rep(c(0.1), 10))
for (i in 1:10) {
  trainSubset <- train[trainSplit != i, 2:ncol(train)]
  testSubset <- train[trainSplit == i, 2:ncol(train)]
  
  model <- naive_bayes(functional_cluster ~ ., data = trainSubset, usekernel = T)
  predictions <- predict(model, testSubset[,-c(ncol(testSubset))])
  table1 <- table(predictions, testSubset$functional_cluster)
  
  classifiedRate = c(classifiedRate, sum(diag(table1))/sum(table1))
}
sprintf("The average correct classification rate: %f", mean(classifiedRate))

# Performance with random guesses
randomClassifiedRate = c()
for (i in 1:100){
  randomPredictions = sample(c(1, 2, 6, 7, 10, 12), nrow(train), replace=T, 
                             prob=c(nrow(train[train$functional_cluster == 1,])/nrow(train), 
                                    nrow(train[train$functional_cluster == 2,])/nrow(train), 
                                    nrow(train[train$functional_cluster == 6,])/nrow(train), 
                                    nrow(train[train$functional_cluster == 7,])/nrow(train),
                                    nrow(train[train$functional_cluster == 10,])/nrow(train),
                                    nrow(train[train$functional_cluster == 12,])/nrow(train)))
  actual = as.vector(train$functional_cluster)
  randomPredictionsTable = table(actual, randomPredictions)
  randomClassifiedRate = c(randomClassifiedRate, sum(diag(randomPredictionsTable))/sum(randomPredictionsTable))
}
sprintf("The average classification rate with random predictions on training data: %f", mean(randomClassifiedRate))

# The performance of our model is not that much better than random predictions...

# Try out model performance with test set
model <- naive_bayes(functional_cluster ~ ., data = train, usekernel = T)
testPredictions <- predict(model, test[,-c(ncol(test))])
functionConfusionMatrix <- table(testPredictions, test$functional_cluster)
testClassifiedRate = sum(diag(functionConfusionMatrix))/sum(functionConfusionMatrix)
sprintf("The correct classification rate on the test data: %f", testClassifiedRate)
