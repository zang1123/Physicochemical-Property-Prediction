#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: MP.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict melting point (MP)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/MP"

# Load MP data

MPdata<-read.table(file.path(outDir, "MP-FP500-MW.txt"), header=T, sep="\t", as.is=T)

dim(MPdata)    # 8648  506  

# 500 fingerprint bits + MW + MP

# There are 6485 training samples and 2163 test samples. Total: 8648

MPdataTraining <-MPdata[1:6485, 1:502]

MPdataTest <-MPdata[6486:8648, 1:502]

dim(MPdataTraining)               # 6485  502

dim(MPdataTest)                      # 2163  502


# Load package e1071 for SVM modeling

library(e1071)

# Use the function svm() to build the SVM model

MPmodel <- svm(MP~., data=MPdataTraining, cost = 9, epsilon = 0.18, gamma = 0.00065)

# Correlation between experimental values and predicted values for training set

MPTrainingCorr<-lm(MPmodel$fitt ~ MPdataTraining$MP)

summary (MPTrainingCorr)

# Predict MP from the test set

MPpred<-predict(MPmodel, MPdataTest)

# Correlation between experimental values and predicted values for test set

MPTestCorr<-lm(MPpred ~ MPdataTest$MP)

summary (MPTestCorr)


setwd("C:/PropertyPrediction/")

save(MPmodel, file="MPmodel.rda")

load("MPmodel.rda")


# Descriptor names for regression modeling

MPfingerprints<-names(MPdata[, c(1:501)])               


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
  modelData<-data[, MPfingerprints]
  #run model
  MPpred<-predict(MPmodel, modelData)
  
  #save output
  write.table(MPpred, file.path(outDir, paste("MPoutput", i, ".txt", sep = "")))
}
