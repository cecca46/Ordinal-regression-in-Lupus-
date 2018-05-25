

require(sva)

#Perform combat on the gene expression data to "remove" batch effect

batch_mitigation <- function(edata, pheno){
    
    btch = pheno$Platform
    #Box plot prior performing quantile norm 
    
    group = pheno$Group
    mod <- model.matrix(~as.factor(group))
    z = ComBat(edata, btch, mod = mod, par.prior=TRUE, prior.plots=TRUE)
    return (z)
}














