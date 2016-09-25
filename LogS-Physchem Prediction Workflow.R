##Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: LogS.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict water solubility (LogS)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/LogS"

# Load LogS data

LogSdata<-read.table(file.path(outDir, "LogS-FP350-MW.txt"), header=T, sep="\t", as.is=T)

dim(LogSdata)    # 2010 356  

# 350 fingerprint bits + MW + LogS

# There are 1507 training chemicals and 503 test chemicals

LogSdataTraining <-LogSdata[1:1507, 1:352]

LogSdataTest <-LogSdata[1508:2010, 1:352]

dim(LogSdataTraining)               # 1507 352

dim(LogSdataTest)                      # 503 352


# Load package e1071 for SVM modeling

library(e1071)

# Use the function svm() to build the SVM model

LogSmodel <- svm(LogS~., data=LogSdataTraining, cost = 260, epsilon = 0.145, gamma = 0.000031)

# Correlation between experimental values and predicted values for training set

LogSTrainingCorr<-lm(LogSmodel$fitt ~ LogSdataTraining$LogS)

summary (LogSTrainingCorr)

# Predict LogS from the test set

LogSpred<-predict(LogSmodel, LogSdataTest)

# Correlation between experimental values and predicted values for test set

LogSTestCorr<-lm(LogSpred ~ LogSdataTest$LogS)

summary (LogSTestCorr)


setwd("C:/PropertyPrediction/")

save(LogSmodel, file="LogSmodel.rda")

load("LogSmodel.rda")


# Descriptor names for regression modeling

LogSfingerprints<-names(LogSdata[, c(1:351)])               


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
  modelData<-data[, LogSfingerprints]
  #run model
  LogSpred<-predict(LogSmodel, modelData)
  
  #save output
  write.table(LogSpred, file.path(outDir, paste("logSoutput", i, ".txt", sep = "")))
}
