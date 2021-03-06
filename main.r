
{
    
# MAIN SCRIPT TO GET DATA, RMA NORMALIZE THE DATA, PERFORM ORDINAL REGRESSION ANALYSIS
# SET UP WORKING ENVIRONMENT 

# install necessary packages if not present yet 
list.of.packages <- c("GEOquery", "affyio","affy", "hgu133a.db","hgu133acdf", "hgu133plus2.db","hgu133plus2cdf","sva","factoextra",
                      "ggplot2","ordinal","cluster","ggfortify","biomaRt","piano", "snowfall", "snow","VennDiagram")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


 
# make sure to set the working directory to the location this function is stored

# setwd("current location")
file_list <- list.files(getwd())
if (!'OrdinalRegression.r' %in% file_list && !'main.r' %in% file_list){
    stop("Please set the working directory to be the one containing this function")
}


# source required functions
source("QuantileNorm.r")
source("PCAplots.r")
source("RMA.r")
source("Combat.r")
source("OrdinalRegression.r")


# set up directories to store results and plots
dir.create("RESULT")
dir.create("PLOT")


# FINISH SETTING UP WORKING ENVIRONMENT -----------------------------



# DOWNLOADING DATASET
# RMA NORMALIZATION
# MERGING TWO DIFFERENT PLATFORMS IN THE DATASET

dataset <- "GSE104948"
get_GEO(dataset)
aggregate_dataset("RESULT/PlatformU133.csv", "RESULT/Platform2Plus.csv", "RESULT/MergedDatasets.csv")
# function that corrects the gene names wrongly changed to DATE by the csv formatter
fix_CSVformat("RESULT/MergedDatasets.csv")

# FINISH DOWNLOADING AND MERGING DATA -------------------------------


# PLOT PCA OF THE MERGED AND NORMALIZED DATASET
edata <- read.csv('RESULT/MergedDatasets.csv')
pheno <- read.csv('MergedDatasetsInfo.csv')
rownames(edata) <- edata[,1]
edata <- edata[,-1]

# saving plot in the RESULT folder
png(filename="PLOT/VariancePlot.png")
variance_plot(edata)
dev.off()
png(filename="PLOT/PCAbyGroup.png")
pca_plot(edata, pheno$Group)
dev.off()
png(filename="PLOT/PCAbyPlatform.png")
pca_plot(edata, pheno$Platform)
dev.off()

# FINISH PLOTTING PCA ------------------------------------------------



# PERFORM QUANTILE NORMALIZATION
# RUN COMBAT TO CORRECT FOR THE PLATFORM TECHNICAL VARIATION CLEAR FROM PREVIOUS PLOTS

png(filename="PLOT/BOXPLOTBEFOREQuantNorm.png")
boxplot( edata, las = 2 )
dev.off()
edata <- t(quantileNormalization(t(edata)))
png(filename="PLOT/BOXPLOTAFTERQuantNorm.png")
boxplot( edata, las = 2 )
dev.off()
# batch = platform
# covariate = condition (Lupus vs Healthy)
batch_corrected <- batch_mitigation(edata, pheno)

# inspect PCA plots again
png(filename="PLOT/PCAbyGroupCombat.png")
pca_plot(batch_corrected, pheno$Group)
dev.off()
png(filename="PLOT/PCAbyPlatformCombat.png")
pca_plot(batch_corrected, pheno$Platform)
dev.off()

# save the batch corrected results for Ordinal regression analysis 
write.csv(batch_corrected, "RESULT/AdjustedLupus.csv")


# FINISH CORRECTING THE DATA ----------------------------------------------


# START ORDINAL REGRESSION ANALYSIS --------------------------------------
# read gene expression matrix and clinical information about patients 

batch_corrected <- read.csv('RESULT/AdjustedLupus.csv')
rownames(batch_corrected) <- batch_corrected[,1]
batch_corrected <- batch_corrected[,-1]
batch_corrected <- t(batch_corrected)
pheno <- read.csv('MergedDatasetsInfo.csv')
# ordinal_result is a list that contains:
# 1 dataframe with the statistics (beta coefficient and pvalues) for all the genes (for GO analysis) and
# 1 dataframe with the statistics (beta coefficient and pvalues) for the significant genes 
# Gene significant is assessed using a default threshold pvalue of 1e-06
# High positive beta means that the gene progressively upregulates with disease progression 
# Small negative beta means that the gene progressively downregulated with disease progression
ordinal_result_list <- ordinal_regression(batch_corrected, pheno) 


all_ordinal_result <- ordinal_result_list$all_ordinal_result
significant_ordinal_result <- ordinal_result_list$significant_ordinal_result
# writing the results of the ordinal analysis: 
# list of significant genes along with pvalue and beta parameter for each gene
write.csv(significant_ordinal_result, "RESULT/SignificantOrdinalRegressionResults.csv")
write.csv(all_ordinal_result, "RESULT/AllOrdinalRegressionResults.csv")
# fix names problem caused by csv formatter
fix_CSVformat("RESULT/SignificantOrdinalRegressionResults.csv")
fix_CSVformat("RESULT/AllOrdinalRegressionResults.csv")



# perform Piano analysis with GO terms on the complete list of results from Ordinal regression
# piano uses 5 methods and pvalues to calculate significance
# a consensus score obteined from all the methods is then computed
# the plots arguments (TRUE by default) saves to the PLOT folder two Venn Diagrams showing the intercepts between the 
# the non directional, directional up, mixed up and non directional, directional down, mixed down respectively
# the method returns a list containing 2 vectors:
# a vector with the significant upregulated bp and 
# a vector with the significant downregulated bp
newList <- run_piano("RESULT/AllOrdinalRegressionResults.csv", plots = TRUE)
write.csv(newList$up,"RESULT/UpRegulatedBP.csv")
write.csv(newList$down,"RESULT/DOWNRegulatedBP.csv")

}

