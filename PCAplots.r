
#Wrapper functions to inspect the data and look for batch effect

#Function to plot the variance explement by the first 10 components
variance_plot <- function(exp_data){
    
    if (class(exp_data[,1]) == 'factor'){
        rownames(exp_data) <- exp_data[, 1]
        exp_data <- exp_data[,-1]
    }
    
    res.pca <- prcomp(exp_data)
    plot(fviz_eig(res.pca))
    summary(res.pca)
    
}


pca_plot <- function(exp_data, phenoCharacter){
    
    if (class(exp_data[,1]) == 'factor'){
        rownames(exp_data) <- exp_data[, 1]
        exp_data <- exp_data[,-1]
    }
    all <- t(exp_data)
    df_pca <- prcomp(all)
    df_out <- as.data.frame(df_pca$x)
    p <-ggplot(df_out,aes(x=PC1,y=PC2,color=phenoCharacter)) + geom_point()
    p <- p + geom_text(aes(label=ifelse(pheno$Samples == 'GSM2810714_H5.Glom.LD_478' | pheno$Samples == 'GSM2810715_H5.Glom.LD_479' | pheno$Samples == 'GSM2810716_H5.Glom.LD_482', as.character(pheno$Samples),"")) ,hjust=0, vjust=0)
    #p <- p +geom_text(aes(label = pheno$Samples),hjust=0, vjust=0)
    plot(p)
    
    
}
pca_plot_2 <- function(exp_data, phenoData){
    
    if (class(exp_data[,1]) == 'factor'){
        rownames(exp_data) <- exp_data[, 1]
        exp_data <- exp_data[,-1]
    }
    prinComp <- cbind(phenoData, pcs$rotation)
    plot(prinComp[, c("stages", "PC1", "PC2", "PC3")], pch = 19, cex = 0.8)
    
    
}

box_plot <- function(exp_data, las = 2){
    
    if (class(exp_data[,1]) == 'factor'){
        rownames(exp_data) <- exp_data[, 1]
        exp_data <- exp_data[,-1]
    }
    boxplot( exp_data, las = las )
}





