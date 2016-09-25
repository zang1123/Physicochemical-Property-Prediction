
#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: BCF.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict bioconcentration factor (BCF)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/BCF"

# Load BCF data

BCFdata<-read.table(file.path(outDir, "BCF-FP200-MW.txt"), header=T, sep="\t", as.is=T)

dim(BCFdata)    # 608    206  

# 200 fingerprint bits + MW + LogBCF

# There are 456 training samples and 152 test samples. Total: 608

BCFdataTraining <-BCFdata[1:456, 1:202]

BCFdataTest <-BCFdata[457:608, 1:202]

dim(BCFdataTraining)               # 456    202

dim(BCFdataTest)                      # 152   202


# Load package e1071 for SVM modeling

library(e1071)

# Use the function svm() to build the SVM model

BCFmodel <- svm(LogBCF~., data=BCFdataTraining, cost = 5300, epsilon=0.10, gamma = 0.000039)

# Correlation between experimental values and predicted values for training set

BCFTrainingCorr<-lm(BCFmodel$fitt ~ BCFdataTraining$LogBCF)

summary (BCFTrainingCorr)

# Predict BCF from the test set

BCFpred<-predict(BCFmodel, BCFdataTest)

# Correlation between experimental values and predicted values for test set

BCFTestCorr<-lm(BCFpred ~ BCFdataTest$LogBCF)

summary (BCFTestCorr)


setwd("C:/PropertyPrediction/")

save(BCFmodel, file="BCFmodel.rda")

load("BCFmodel.rda")


# Descriptor names for regression modeling

BCFfingerprints<-names(BCFdata[, c(1:201)])               

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
  modelData<-data[, BCFfingerprints]
  #run model
  BCFpred<-predict(BCFmodel, modelData)
  
  #save output
  write.table(BCFpred, file.path(outDir, paste("BCFoutput", i, ".txt", sep = "")), sep = "\t")
}
