#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: LogP.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict octonal-water partition coefficient (LogP)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/LogP"

# Load LogP data

LogPdata<-read.table(file.path(outDir, "LogP-FP600-MW.txt"), header=T, sep="\t", as.is=T)

dim(LogPdata)    # 14193   606  

# 600 fingerprint bits + MW + LogP

# There are 11371 training samples and 2837 test samples. Total: 14208

LogPdataTraining <-LogPdata[1:11360, 1:602]

LogPdataTest <-LogPdata[11361:14193, 1:602]

dim(LogPdataTraining)               # 11360      602

dim(LogPdataTest)                      # 2833         602


# Load support vector machine (SVM) package

library(e1071)


# Use the function svm() to build the SVM model

LogPmodel <- svm(LogP~., data=LogPdataTraining, cost = 150, epsilon = 0.025, gamma = 0.00014)

# Correlation between experimental values and predicted values for training set

LogPTrainingCorr<-lm(LogPmodel$fitt ~ LogPdataTraining$LogP)

summary (LogPTrainingCorr)


# Predict LogS from the test set

LogPpred<-predict(LogPmodel, LogPdataTest)


# Correlation between experimental values and predicted values for test set

LogPTestCorr<-lm(LogPpred ~ LogPdataTest$LogP)

summary (LogPTestCorr)


setwd("C:/PropertyPrediction/")

save(LogPmodel, file="LogPmodel.rda")

load("LogPmodel.rda")


# Descriptor names for regression modeling

LogPfingerprints<-names(LogPdata[, c(1:601)])               


filenames<-c("DSSTox-QSAR1-FP-MW.txt", 
             "DSSTox-QSAR2-FP-MW.txt",
             "DSSTox-QSAR3-FP-MW.txt",
             "DSSTox-QSAR4-FP-MW.txt",
             "DSSTox-QSAR5-FP-MW.txt",
             "DSSTox-QSAR6-FP-MW.txt",
             "DSSTox-QSAR7-FP-MW.txt",
             "DSSTox-QSAR8-FP-MW.txt",
             "DSSTox-QSAR9-FP-MW.txt",
             "DSSTox-QSAR10-FP-MW.txt",
             "DSSTox-QSAR11-FP-MW.txt",
             "DSSTox-QSAR12-FP-MW.txt",
             "DSSTox-QSAR13-FP-MW.txt",
             "DSSTox-QSAR14-FP-MW.txt",
             "DSSTox-QSAR15-FP-MW.txt")


#define input direcory
inDir<-"C:/PropertyPrediction/Descriptors"

for(i in 1:length(filenames)){
  data<-read.table(file.path(inDir, filenames[i]), header=T, sep="\t", as.is=T)
  #select appropriate fingerprints
  modelData<-data[, LogPfingerprints]
  #run model
  LogPpred<-predict(LogPmodel, modelData)
  
  #save output
  write.table(LogPpred, file.path(outDir, paste("logPoutput", i, ".txt", sep = "")))
}
