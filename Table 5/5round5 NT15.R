label = read.csv("Class angio and negative.csv", header = TRUE) 
data2 = data.frame(PCP,Class = label)

Pos = subset(data2, Class == 'Antiangio')
Neg = subset(data2, Class == 'Negative')
nPos = nrow(Pos)
nNeg = nrow(Neg)

m= 10
ACCtr  <- matrix(nrow = m, ncol = 1)
SENStr  <- matrix(nrow = m, ncol = 1)
SPECtr  <- matrix(nrow = m, ncol = 1)
MCCtr  <- matrix(nrow = m, ncol = 1)
AUCtr  <- matrix(nrow = m, ncol = 1)
ACCts  <- matrix(nrow = m, ncol = 1)
SENSts  <- matrix(nrow = m, ncol = 1)
SPECts  <- matrix(nrow = m, ncol = 1)
MCCts  <- matrix(nrow = m, ncol = 1)
AUCts  <- matrix(nrow = m, ncol = 1)

for (i in 1:m){
#######  Dividing Training and Testing sets on positive and negative classes
sample1 <- c(sample(1:99,80))
sample2 <- c(sample(1:101,80))
  train1  <- Pos[sample1,] ####Positive set for training
  train2  <- Neg[sample2,] ####Negative set for training
  test1 <-   Pos[-sample1,]    ####Positive set for testing
  test2 <-   Neg[-sample2,]    ####Negative set for testing 
  internal <- rbind(train1,train2) ####combining for internal set
  external <- rbind(test1,test2)    ####combining for external set

######### External validation over 5CV
tunegrid <- expand.grid(.mtry=c(1:10), .ntree=seq(100,500,100))
RFmodel <- train(Class~., data=internal , method=customRF, metric=c("Accuracy"), tuneGrid=tunegrid, trControl=trainControl(method= "repeatedcv", number=5,repeats=5))
Resultcv = RFmodel$ finalModel$ confusion [,1:2]
pred=prediction(RFmodel$ finalModel$ votes[,2],internal[,ncol(internal)])
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUCtr[i,] = perf_AUC@y.values[[1]]

################### External validation
predcv = table(external$Class, predict(RFmodel, external))  ###### Prediction on external set
Resultext <- rbind(predcv[1], predcv[3],predcv[2], predcv[4]) ###### Reporting TN,FP,FN,TP
pred= prediction(predict(RFmodel ,external,type = "prob")[,2],external[,ncol(external)])
perf_AUC=performance(pred,"auc") #Calculate the AUC value
AUCts[i,] = perf_AUC@y.values[[1]]
################### Performance report
data = Resultcv
	ACCtr[i,] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENStr[i,]  =  (data[1]/(data[1]+data[2]))*100
	SPECtr[i,] = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[3]*data[2])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCtr[i,]  = MCC1/MCC4
data = Resultext
	ACCts[i,] = (data[1]+data[4])/(data[1]+data[2]+data[3]+data[4])*100
	SENSts[i,]  =  (data[1]/(data[1]+data[2]))*100
	SPECts[i,] = (data[4])/(data[3]+data[4])*100
	MCC1      = (data[1]*data[4]) - (data[3]*data[2])
	MCC2      =  (data[4]+data[2])*(data[4]+data[3])
	MCC3      =  (data[1]+data[2])*(data[1]+data[3])
	MCC4	=  sqrt(MCC2)*sqrt(MCC3)
	MCCts[i,]  = MCC1/MCC4
}

result = data.frame (ACCtr,SPECtr,SENStr,MCCtr,AUCtr,ACCts,SPECts,SENSts,MCCts,AUCts)
result = na.omit(result)

Mean  <- matrix(nrow = 10, ncol = 1)
SD  <- matrix(nrow = 10, ncol = 1)
for (i in 1:10){
Mean[i,] = mean(result[,i])
SD[i,] =  sd(result[,i])
}
round(data.frame (Mean,SD),2)

