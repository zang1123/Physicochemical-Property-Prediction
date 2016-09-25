#Author: Dan Zang, PhD
#Affiliation: Integrated Laboratory Systems, Inc. 
#Contact: dzang@ils-inc.com
#Title: BP.R 
#Purpose: this code takes an input of PaDEL descriptors and uses the SVM model from e1071 #package to predict boiling point (BP)
#Notes: genetic algorithm was used to determine the most relevant molecular fingerprints from #the PaDEL descriptors. 

#define output direcory

outDir<-"C:/PropertyPrediction/BP"

# Load BP data

BPdata<-read.table(file.path(outDir, "BP-FP400-MW.txt"), header=T, sep="\t", as.is=T)

dim(BPdata)    # 5432  406  

# 400 fingerprint bits + MW + BP

# There are 4074 training samples and 1358 test samples. Total: 5432

BPdataTraining <-BPdata[1:4074, 1:402]

BPdataTest <-BPdata[4075:5432, 1:402]

dim(BPdataTraining)               # 4074 402

dim(BPdataTest)                      # 1358 402


# Load package e1071 for SVM modeling

library(e1071)

# Use the function svm() to build the SVM model

BPmodel <- svm(BP~., data=BPdataTraining, cost = 9, epsilon=0.012, gamma = 0.0010)

# Correlation between experimental values and predicted values for training set

BPTrainingCorr<-lm(BPmodel$fitt ~ BPdataTraining$BP)

summary (BPTrainingCorr)

# Predict BP from the test set

BPpred<-predict(BPmodel, BPdataTest)

# Correlation between experimental values and predicted values for test set

BPTestCorr<-lm(BPpred ~ BPdataTest$BP)

summary (BPTestCorr)


setwd("C:/PropertyPrediction/")

save(BPmodel, file="BPmodel.rda")

load("BPmodel.rda")


# Descriptor names for regression modeling

BPfingerprints<-names(BPdata[, c(1:401)])               


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
  modelData<-data[, BPfingerprints]
  #run model
  BPpred<-predict(BPmodel, modelData)
  
  #save output
  write.table(BPpred, file.path(outDir, paste("BPoutput", i, ".txt", sep = "")))
}
