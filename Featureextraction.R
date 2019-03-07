#######set directory
setwd('D:\\Peptide prediction\\Anti-Angiogenic Peptides\\backup')
#######Load package
library(randomForest)
library(protr)
library(seqinr)
library(C50)
library(RWeka)
library(Interpol)
library(caret)
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

####### Angio vs Negative
x <- read.fasta('antiangio and negative.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
PCP  <- matrix(nrow = m, ncol = 531)
pse = 1
weightpse = 0.9
ampse = 1
weightampse = 0.9

PAAC <- matrix(nrow = m, ncol = 20 + pse)
AAC <- t(sapply(A, extractAAC))
DPC <- t(sapply(A, extractDC))

for(i in 1:m){ 
x = A[[i]][1]
b = AAdescriptor(x, 531,2)
PCP[i,] = Interpol(b, 531,"linear")
}

for(i in 1:m){ 
PAAC[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = weightpse, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

col = 20+ 2*ampse
APAAC  <- matrix(nrow = length(A), ncol = col)
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = ampse, w = weightampse, customprops = NULL)
}

####### Angio vs Negative NTCT
x <- read.fasta('angio NT15.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative NTCT.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
PCP  <- matrix(nrow = m, ncol = 531)
pse = 2
weight = 0.1
PAAC <- matrix(nrow = m, ncol = 20 + pse)

AAC <- t(sapply(A, extractAAC))
DPC <- t(sapply(A, extractDC))
CTDC <- t(sapply(A, extractCTDC))
CTDT <- t(sapply(A, extractCTDT))
CTDD <- t(sapply(A, extractCTDD))

for(i in 1:m){ 
x = A[[i]][1]
b = AAdescriptor(x, 531,2)
PCP[i,] = Interpol(b, 531,"linear")
}

for(i in 1:m){ 
PAAC[i, ] = extractPAAC(A[[i]][1],lambda = pse, w = weight, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"))
}

col = 20+ 2*3
APAAC  <- matrix(nrow = length(A), ncol = col)
for (i in 1:length(A)){
APAAC[i,] = extractAPAAC(A[[i]][1],lambda = 3, w = 0.2, customprops = NULL)
}
