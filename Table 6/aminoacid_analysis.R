
library(RWeka)
library(caret)
library(randomForest)
library(Interpol)
library(protr)
library(seqinr)
library(Peptides)

#######Extract feature
x <- read.fasta('antiangio and negative.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
AAC <- t(sapply(A, extractAAC))
data = data.frame(AAC, Class = label)

meanX  <- matrix(nrow = ncol(data), ncol = 1)
sdX  <- matrix(nrow = ncol(data), ncol = 1)
meanY  <- matrix(nrow = ncol(data), ncol = 1)
sdY  <- matrix(nrow = ncol(data), ncol = 1)
p.value  <- matrix(nrow = ncol(data), ncol = 1)
m = ncol(data)
k = m-1

for(i in 1:k){ 
	X <- subset(data.frame(data[,i],Activities =data[,m]), Activities == 'Antiangio')
	Y <- subset(data.frame(data[,i],Activities =data[,m]), Activities == 'Negative')
	meanX[i,]  = mean(as.numeric(X[,1]))
	sdX[i,]  = sd(as.numeric(X[,1]))
	meanY[i,]  = mean(as.numeric(Y[,1]))
	sdY[i,]  = sd(as.numeric(Y[,1]))
	p.value[i,]  = t.test(X[,1],Y[,1])$ p.value
}

stat = data.frame (meanX,sdX,meanY,sdY,p.value)

##########MDGI calculation

internal = data

ind= c(2,3,5,7,9,11,13,15,17,20)
n = ncol(internal)-1
gini = matrix(nrow = n, ncol = 10)
meangini = matrix(nrow = n, ncol = 1)

for (i in 1:10){
RF<-randomForest(Class~.,data=internal,ntree=100,mtry=ind[i],importance=TRUE)
gini[,i] = RF$ importance[,4]
}

for (i in 1:n){
meangini[i,] = mean(gini[i,])
}
