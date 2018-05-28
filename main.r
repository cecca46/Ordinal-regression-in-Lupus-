

#MAIN SCRIPT TO GET DATA, RMA NORMALIZE THE DATA, PERFORM ORDINAL REGRESSION ANALYSIS


#SET UP WORKING ENVIRONMENT 

#Save current user directory
USER_WD <- getwd()

#Change Working directory to the one containing scripts necessary to perform the analysis
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
#Source required functions
source("QuantileNorm.r")
source("PCAplots.r")
source("RMA.r")
source("Combat.r")
source("OrdinalRegression.r")

#Set up directory to store results
dir.create("RESULT")

#FINISH SETTING UP WORKING ENVIRONMENT -----------------------------




#DOWNLOADING DATASET
#RMA NORMALIZATION
#MERGING TWO DIFFERENT PLATFORMS IN THE DATASET

dataset <- "GSE104948"
get_GEO(dataset)
aggregate_dataset("RESULT/PlatformU133.csv", "RESULT/Platform2Plus.csv", "RESULT/MergedDatasets.csv")

#FINISH DOWNLOADING AND MERGING DATA -------------------------------


#PLOT PCA OF THE MERGED AND NORMALIZED DATASET
edata <- read.csv('RESULT/MergedDatasets.csv')
pheno <- read.csv('MergedDatasetsInfo.csv')
rownames(edata) <- edata[,1]
edata <- edata[,-1]

#Saving plot in the RESULT folder
pdf(file="RESULT/VariancePlot.pdf")
variance_plot(edata)
dev.off()
pdf(file="RESULT/PCAbyGroup.pdf")
pca_plot(edata, pheno$Group)
dev.off()
pdf(file="RESULT/PCAbyPlatform.pdf")
pca_plot(edata, pheno$Platform)
dev.off()

#FINISH PLOTTING PCA ------------------------------------------------



#PERFORM QUANTILE NORMALIZATION
#RUN COMBAT TO CORRECT FOR THE PLATFORM TECHNICAL VARIATION CLEAR FROM PREVIOUS PLOTS

pdf(file="RESULT/BOXPLOTBEFOREQuantNorm.pdf")
boxplot( edata, las = 2 )
dev.off()
edata <- t(quantileNormalization(t(edata)))
pdf(file="RESULT/BOXPLOTAFTERQuantNorm.pdf")
boxplot( edata, las = 2 )
dev.off()
#Batch = platform
#Covariate = condition (Lupus - Healthy)
batch_corrected <- batch_mitigation(edata, pheno)

#Inspect PCA plots again
pdf(file="RESULT/PCAbyGroupCombat.pdf")
pca_plot(batch_corrected, pheno$Group)
dev.off()
pdf(file="RESULT/PCAbyPlatformCombat.pdf")
pca_plot(batch_corrected, pheno$Platform)
dev.off()

#Save the batch corrected results for Ordinal regression analysis 
write.csv(batch_corrected, "RESULT/AdjustedLupus.csv")


#FINISH CORRECTING THE DATA ----------------------------------------------


#START ORDINAL REGRESSION ANALYSIS --------------------------------------
#Read gene expression matrix and clinical information about patients 

batch_corrected <- read.csv('RESULT/AdjustedLupus.csv')
rownames(batch_corrected) <- batch_corrected[,1]
batch_corrected <- batch_corrected[,-1]
batch_corrected <- t(batch_corrected)
pheno <- read.csv('MergedDatasetsInfo.csv')
# ordinal_result is a dataframe that contains a number of significant genes along with their beta coefficient 
# the genes are ordered according to the value of the beta coefficient
# high positive beta means that the gene progressively upregulates with disease progression 
# small negative beta means that the gene progressively downregulated with disease progression
ordinal_result <- ordinal_regression(batch_corrected, pheno) 
write.csv(ordinal_result, "RESULT/OrdinalRegressionResults.csv")








#AT THE END! RESET USER WORKING ENVIRONMENT AS BEFORE RUNNING THIS FUNCTION 
setwd(USER_WD)
