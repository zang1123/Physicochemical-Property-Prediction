#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: VP.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict vapor pressure (VP)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/VP"

# Load MP data

VPdata<-read.table(file.path(outDir, "VP-FP350-MW.txt"), header=T, sep="\t", as.is=T)

dim(VPdata)    # 2713    356  

# 350 fingerprint bits + MW + LogVP

# There are 2034 training samples and 679 test samples. Total: 2713

VPdataTraining <-VPdata[1:2034, 1:352]

VPdataTest <-VPdata[2035:2713, 1:352]

dim(VPdataTraining)               # 2034 352

dim(VPdataTest)                      # 679 352


# Load package e1071 for SVM modeling

library(e1071)

# Use the function svm() to build the SVM model

VPmodel <- svm(LogVP~., data=VPdataTraining, cost = 170, epsilon = 0.07, gamma = 0.00011)

# Correlation between experimental values and predicted values for training set

VPTrainingCorr<-lm(VPmodel$fitt ~ VPdataTraining$LogVP)

summary (VPTrainingCorr)

# Predict VP from the test set

VPpred<-predict(VPmodel, VPdataTest)

# Correlation between experimental values and predicted values for test set

VPTestCorr<-lm(VPpred ~ VPdataTest$LogVP)

summary (VPTestCorr)


setwd("C:/PropertyPrediction/")

save(VPmodel, file="VPmodel.rda")

load("VPmodel.rda")


# Descriptor names for regression modeling

VPfingerprints<-names(VPdata[, c(1:351)])               


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
  modelData<-data[, VPfingerprints]
  #run model
  VPpred<-predict(VPmodel, modelData)
  
  #save output
  write.table(VPpred, file.path(outDir, paste("VPoutput", i, ".txt", sep = "")))
}
