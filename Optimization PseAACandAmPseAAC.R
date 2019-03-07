#######set directory
setwd('D:\\Peptide prediction\\Anti-Angiogenic Peptides\\backup')
#######Load package
library(Interpol)
library(caret)
library(randomForest)
library(cvTools) 
library(Metrics)
library(MASS)
library(pls)
library(Interpol)
library(protr)
library(seqinr)
library(Peptides)
library(AUC)
library(ROCR)

customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes

#######Optimizing PseAAC
x <- read.fasta('angio CT15.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative NTCT.csv", header = TRUE) 

lam = 10
w = 10
A <- x[(sapply(x, protcheck))]
m = length(A)
index = seq(0.1,1, 0.1)
ACC  <- matrix(nrow = length(index), ncol = lam)

for(pse in 1:lam){
paac <- matrix(nrow = m, ncol = 20 + pse)
for(weight in 1:w){
for(i in 1:m){ 
paac[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = index[weight], customprops = NULL)
}
internal= data.frame (paac, Class = label)
M = train(Class ~ ., data = internal, 
          method = "gaussprRadial",
           trControl = trainControl(method = "cv", number = 10, savePredictions = TRUE),
           preProc = c("center", "scale"))

data = table(M$pred[,1],M$pred[,2])
ACC[weight,pse] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
}}


#######Optimizing Am-PseAAC
x <- read.fasta('angio CT15.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative NTCT.csv", header = TRUE) 

lam = 10
w = 10
A <- x[(sapply(x, protcheck))]
m = length(A)
index = seq(0.1,1, 0.1)
ACC  <- matrix(nrow = length(index), ncol = lam)

for(pse in 1:lam){
col = 20+ 2*pse
APAAC  <- matrix(nrow = length(A), ncol = col)
for(weight in 1:w){
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = pse, w = index[weight], customprops = NULL)
}
internal= data.frame (APAAC, Class = label)
M = train(Class ~ ., data = internal, 
          method = "gaussprRadial",
           trControl = trainControl(method = "cv", number = 10, savePredictions = TRUE),
           preProc = c("center", "scale"))

data = table(M$pred[,1],M$pred[,2])
ACC[weight,pse] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
}}

