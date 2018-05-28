
#Function to perform ordinal regression.
#input: gene expression matrix in the form of dataframe with samples in rows and genes in columns, dataframe or vector containing the ordinal variable of interest 
#output: list of upregulated and downregulated PAGs (phenotype associated genes)

ordinal_regression <- function(exp, pheno){
    #read from the clinical information the ordinal categories (stages)
    #make the ordinal category a column of the gene matrix
    ordinal_variable <- pheno$stages
    exp <- cbind(exp, ordinal_variable)
    #adjust the columns name 
    colnames(exp) <- gsub(pattern = "-",replacement = ".",x = colnames(exp))
    colnames(exp) <- gsub('\"', "", colnames(exp), fixed = TRUE)
    colnames(exp) <- gsub(pattern = "@", replacement = "",x = colnames(exp))
    
    exp <- as.data.frame(exp)
    #all columns (features) of exp excluding the response variable
    varNames <- colnames(exp)
    varNames <- varNames[!varNames %in% c("ordinal_variable")]
    varNames <- gsub(pattern = "-",replacement = ".",x = varNames)
    
    #make the output variable as numeric
    #group the conditions of the samples to increase their number per individual ordinal category
    #Ordinal category 1: Controls
    #Ordinal category 2: Stage 1 and Stage 2
    #Ordinal category 3: Stage 3A and Stage 3B
    #Ordinal category 4: Stage 4 and Stage 5
    exp$ordinal_variable <- as.numeric(exp$ordinal_variable)
    exp$ordinal_variable[exp$ordinal_variable == '2' | exp$ordinal_variable == '3'] <- '2'
    exp$ordinal_variable[exp$ordinal_variable == '4' | exp$ordinal_variable == '5'] <- '3'
    exp$ordinal_variable[exp$ordinal_variable == '6' | exp$ordinal_variable == '7'] <- '4'
    table(exp$ordinal_variable)
    #make response variable an ordered factor
    exp$ordinal_variable <- as.ordered(exp$ordinal_variable)
    class(exp$ordinal_variable)
    #create variables to keep truck of the regression results
    wald <- list()
    genes <- NULL
    betas <- NULL
    
    #fit a regression model for each gene
    for (i in 1:length(varNames)){
        #print (i)
        predictor <- varNames[i]
        form <- as.formula(paste("ordinal_variable~", paste(gsub('\"', "", predictor, fixed = TRUE)),collapse = ""))
        model1 <- clm(form, data = exp, link = 'probit')
        (summ <- summary(model1))
        #save in 3 ordered list (indexed by the same index for the same sample) the gene name, the significance value (p-value from Wald statistics) and beta coefficients
        genes <- c(genes, predictor)
        wald[[predictor]] <- coef(summ)[1:4,4]
        betas <- c(betas, summ$beta)
    }
    #check the significance of each gene
    #a gene is called significant if all the fitted values for that gene are <= threshold
    threshold <- 1e-06 #starting threshold value
    while (1){ #repeat this loop until exit (break) condition is met 
        significant <- find_significant(wald, threshold)
        #if there are too few significant genes, increase the threshold by factor of 10
        if (length(significant) >= 0 & length(significant) <= 400)
            threshold <- threshold * 1e1
        #if there are too many significant results, decrease the threashold by factor of 10
        else if (length(significant) >= 1000){
            threshold <- threshold * 1e-1
        } 
        #if number of significant results is between 400 and 1000, exit the loop
        else
            break()
    }
    #create dataframe for results
    ordinal_result <- data.frame(genes, betas)
    ordinal_result <- ordinal_result[rownames(ordinal_result) %in% significant, ]
    df <- ordinal_result[order(-ordinal_result$betas),]
    return (df)
    #upPAG <- ordinal_result[ordinal_result$betas>0,]
    #downPAG <- ordinal_result[ordinal_result$betas<0,]
    
}
#Function to identify the significant genes according to the dynamic value of theshold
find_significant <- function(pval, threshold){
    
    significant <- NULL
    for(i in 1:length(pval)){
        keep <- TRUE
        for (j in 1:length(pval[[i]])){
            if (pval[[i]][j] >= threshold){
                #if any of the fitted parameters is greated or equal than the threshold than the gene is not considered significant
                #hence not to keep
                keep <- FALSE
            }
        }
        if (keep == TRUE){
            #insert gene name in the list of significant genes 
            significant <- c(significant, names(pval[[i]])[4])
        }
    }
    #list of significant genes 
    return (significant)
    
}
