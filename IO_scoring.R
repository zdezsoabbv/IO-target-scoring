library(randomForest)
library(pROC)
library(caret)
library(foreach)

load("obj.Rdata")

###### cross-validation for 10,000 models by randomly sampling negative set 
Nrun=10000
set.seed(5)
auc=NULL
sens=rep(0,78)
spec=rep(0,78)
for(nr in 1:Nrun)
{ #print(nr)
  neg_train<-sample(neg, length(pos_train), replace = FALSE)
 # neg_train=tneg
  trainData <- features[c(pos_train,neg_train),]
  trainAnnot <- as.factor(c(rep("pos",length(pos_train)),rep("neg",length(neg_train))))  
  ctrl<-trainControl(method="CV",classProbs = TRUE,number=10, savePredictions = T,summaryFunction = twoClassSummary,returnResamp="all")
  tr<-train(trainData,trainAnnot,method="rf",ntree=500,metric="ROC",trControl=ctrl,importance=T)
  selectedIndices <- tr$pred$mtry == tr$finalModel$mtry
  pr=as.character(tr$pred$pos[selectedIndices])
  pr[pr=="neg"]=0
  pr[pr=="pos"]=1
  pr=as.numeric(pr)
  a=roc(tr$pred$obs[selectedIndices],pr)
  thrs <- coords(a, c(-Inf, sort(ctfs), Inf), input="threshold", ret=c("specificity", "sensitivity"))
  auc=c(auc,as.numeric(a$auc))
  spec=spec+thrs$specificity
  sens=sens+thrs$sensitivity
}
spec=spec/Nrun
sens=sens/Nrun

##### cross-validation of model based on negative set

neg_train=tneg
trainData <- features[c(pos_train,neg_train),]
trainAnnot <- as.factor(c(rep("pos",length(pos_train)),rep("neg",length(neg_train))))  
ctrl<-trainControl(method="CV",classProbs = TRUE,number=10, savePredictions = T,summaryFunction = twoClassSummary,returnResamp="all")
tr<-train(trainData,trainAnnot,method="rf",ntree=500,metric="ROC",trControl=ctrl,importance=T)
selectedIndices <- tr$pred$mtry == tr$finalModel$mtry
pr=as.character(tr$pred$pos[selectedIndices])
pr[pr=="neg"]=0
pr[pr=="pos"]=1
pr=as.numeric(pr)
a=roc(tr$pred$obs[selectedIndices],pr)


####### prediction of all genes 

model=vector("list",Nrun)
pred=vector("list",Nrun)
predictData <- features

foreach(i=1:Nrun) %do% 
  { neg_train<-sample(neg, length(pos_train), replace = FALSE)
    trainData <- features[c(pos_train,neg_train),]
    trainAnnot <- as.factor(c(rep(1,length(pos_train)),rep(0,length(neg_train))))
    bestMtry=tuneRF(trainData, trainAnnot, ntreeTry=1000,stepFactor=1.2, improve=0.01, trace=FALSE, plot=FALSE, doBest=FALSE) 
    Mtry=bestMtry[min(which(bestMtry[,2]==min(bestMtry[,2]))),1]
    RF<-randomForest(trainData,trainAnnot,importance=T,ntree=1000,mtry=Mtry)
    model[[i]]<-RF
    pred[[i]]<-predict(RF,predictData,type="prob")
  }

final_prob <- 0 
final_imp <-0 
for(i in 1:length(pred)) { final_prob<- final_prob+pred[[i]] 
final_imp<-final_imp+model[[i]]$importance[,4]
} 
final_prob<- final_prob/length(pred) 
final_imp<-final_imp/length(pred)



