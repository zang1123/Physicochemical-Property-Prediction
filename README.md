# Physicochemical-Property-Prediction

Models for predicting six physicochemical properties for environmental chemicals. The properties are: octonol-water partition coefficient (LogP), water solubility (LogS), melting point (MP), boiling point (BP), vapor pressure (LogVP), and bioconcentration factor (BCF).

Purpose: this code takes an input of binary values from PaDEL molecular fingerprints and uses the model of support vector machine (SVM) to predict the property. The Physchem Prediction Workflow.R is the R code for generating the model.  The modelAR.rda is the SVM model which can be reused to predict the property. BCF-FP200-Mw.txt, BP-FP400-MW.txt, LogP-FP600-MW.txt, LogS-FP350-MW, MP-FP500-MW and VP-FP350-MW.txt contain training and test data for the six properties, respectively.​

This project has been funded in whole or in part with Federal funds from the National Institute of Environmental Health Sciences, National Institutes of Health, Department of Health and Human Services, under Contract No. HHSN273201500010C.

