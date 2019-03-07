#######Load package
library(Interpol)
library(protr)
library(seqinr)

#######Extract AAC DPC TPC
x <- read.fasta('antiangio and negative.fasta', seqtype="AA", as.string = TRUE)
label = read.csv("Class angio and negative.csv", header = TRUE) 
A <- x[(sapply(x, protcheck))]
m = length(A)
PCP  <- matrix(nrow = m, ncol = 531)

for(i in 1:m){ 
x = A[[i]][1]
b = AAdescriptor(x, 531,2)
PCP[i,] = Interpol(b, 531,"linear")
}

##########MDGI calculation

internal = data.frame(PCP, Class = label)

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
