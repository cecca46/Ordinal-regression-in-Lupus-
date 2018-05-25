
#Download GEO dataset 
#Perfom RMA normalization separately on the 2 different platforms 
#Merge the normalized result into a single expression matrix

rm(list = ls())
library(GEOquery)
library(affyio)
library(affy)
library(hgu133a.db)
library(hgu133acdf)
library(hgu133plus2.db)
library(hgu133plus2cdf)
library(sva)
require(factoextra)
require(ggplot2)



#Getting the GEO dataset and extracting sample expression
get_GEO <- function(geoFile){
    
    main_path <- "~/Desktop/"
    # Set working directory for download
    setwd(main_path)
    getGEOSuppFiles(geoFile)
    setwd(paste0(main_path,geoFile,collapse = "/"))
    untar("GSE104948_RAW.tar", exdir = "data")
    cels = list.files("data/", pattern = "gz")
    # sometiles, it is 'CEL', you need to check it first
    setwd("~/Desktop/GSE104948/data/")
    f <- list.files(pattern = "CEL.gz")
    ff <- split(f, sapply(f, function(x) read.celfile.header(x)$cdfName))
    chip_plus2 <- ff[1]
    chip_u133A <- ff[2]
    chip_plus2 <- as.character(unlist(chip_plus2))
    chip_u133A <- as.character(unlist(chip_u133A))
    setwd(paste0(main_path,geoFile,collapse = "/"))
    sapply(paste("data", cels, sep = "/"), gunzip)
    chip_plus2 <- gsub(".gz", "", chip_plus2)
    chip_u133A <- gsub(".gz", "", chip_u133A)
    #cels = list.files("data/", pattern = "CEL")
    chip_plus2 <- chip_plus2[grep('LD|SLE', chip_plus2)]
    chip_u133A <- chip_u133A[grep('LD|SLE', chip_u133A)]
    #perform RMA norm on platform 133 and save results to csv
    RMA_norm(chip_u133A, "PlatformU133.csv", hgu133a.db)
    #perform RMA norm on platform plus2 and save results to csv
    RMA_norm(chip_plus2, "Platform2Plus.csv", hgu133plus2.db)
}


#Perform RMA norm and saves the results in two .csv files according to the platform used to process the samples
RMA_norm <- function(samples, csvName, database){
    
    setwd("~/Desktop/GSE104948/data/")
    raw.data = ReadAffy(verbose = FALSE, filenames = samples)
    # perform RMA normalization (log2)
    data.rma.norm = rma(raw.data)
    rma = exprs(data.rma.norm)
    # Take a look at the result (first 5 rows and 5 columes)
    rma[1:5, 1:5]
    tt = cbind(row.names(rma), rma)
    colnames(tt) = c("ProbID", sub(".cel", "", colnames(rma), ignore.case = TRUE))
    rownames(tt) = NULL
    tt[1:5, 1:5]
    #Mapping gene symbol to entrezId
    annot <- select(x = database, 
                    keys    = tt[,1] ,
                    columns = c("ENTREZID", 'SYMBOL'))
    
    comb = merge(annot, tt, by.x = "PROBEID", by.y = "ProbID")
    comb[1:5, 1:5]
    
    names <- unique(comb$SYMBOL)
    ids <- unique(comb$ENTREZID)
    names <- as.character(names)
    ids <- as.numeric(ids)
    geneNames <- data.frame(names, ids)
    colnames(geneNames) = c("geneName", "entryID")
    
    # If multiple probe sets corresponded to the same gene, then the expression
    # values of these probe sets were averaged.
    
    comb2 <- subset(comb, select = -c(PROBEID, SYMBOL))
    comb2 <- data.frame(lapply(comb2, as.character), stringsAsFactors = FALSE)
    comb2 <- data.frame(lapply(comb2, as.numeric), stringsAsFactors = FALSE)
    out <- aggregate(. ~ ENTREZID, data = comb2, mean)
    out = format(out, digits = 5)
    out[1:5, 1:5]
    out$ENTREZID = as.numeric(out$ENTREZID)
    out2 = merge(geneNames, out , by.x = "entryID", by.y = "ENTREZID")
    out2[1:5, 1:5]
    #Writing result for one of the platform to csv
    write.csv(out2, paste0("~/Desktop/",csvName, collapse = ""))
    
    
}

#Function that reads the two normalised datasets and merge them together
aggregate_dataset <- function(pathToData1, pathToData2, pathToResult){
    
    x <- read.csv(pathToData1)
    y <- read.csv(pathToData2)
    z <- merge(x, y, by="geneName", all.x = T)
    z <- aggregate(. ~ geneName, data = z, mean)
    write.csv(z, pathToResult)
    
    #reding clinical info and subselecting the samples for which clinical info are available
    clin <- read.csv('~/Desktop/MergedDatasetsInfo.csv')
    samples <- as.character(clin$Samples)
    rownames(z) <- z[,1]
    z <- z[,-1]
    edata <- z[,names(z) %in% samples]
    write.csv(edata, pathToResult)
    
}
    

#In line code to start analysis 

get_GEO("GSE104948")
aggregate_dataset("~/Desktop/PlatformU133.csv", "~/Desktop/Platform2Plus.csv", "~/Desktop/MergedDatasets.csv")



