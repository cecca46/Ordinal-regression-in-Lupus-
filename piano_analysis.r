
#function to run piano 
#arguments: gene statistics (pvalues, tvalues and FC to use according to the methods in a dataframe), gene set collection, statistical methods to run and adjustment method.
#default adjastment:FDR


#Piano analysis with GO terms
library(biomaRt)
library(piano)
library(snowfall)
library(snow)
library(VennDiagram)



run_piano <- function(path, plots = TRUE){
    
    
    stat <- read.csv(path)
    rownames(stat) <- stat$X
    stat <- stat[,-1]
    stat <- stat[,-1]
    names(stat) <- c("beta","pval")
    pstat <- subset(stat, select = pval)
    fold_change <- subset(stat, select = beta)
    
    
    # Select ensembl database and hsapiens_gene_ensembl:
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl" ,host="www.ensembl.org")
    # Map human microarray gene names IDs to GO:
    mapGO <- getBM(attributes=c("external_gene_name", "name_1006"), mart = ensembl)
    # Remove blanks ("")
    mapGO <- mapGO[mapGO[,2]!="",]
    mapGO[1:10,]
    #load GSC to run piano
    geneSetCollections <- loadGSC(mapGO)
    
    #run piano different methods
    gsaRes1 <- runGSA(pstat, geneSetStat="fisher", directions = fold_change, gsc=geneSetCollections, adjMethod = 'fdr',nPerm = 1000, ncpus = 4)
    gsaRes2 <- runGSA(pstat, geneSetStat="stouffer", directions = fold_change, gsc=geneSetCollections, adjMethod = 'fdr',nPerm = 1000, ncpus = 4)
    gsaRes3 <- runGSA(pstat, geneSetStat="reporter", directions = fold_change, gsc=geneSetCollections, adjMethod = 'fdr',nPerm = 1000, ncpus = 4)
    gsaRes4 <- runGSA(pstat, geneSetStat="tailStrength", directions = fold_change, gsc=geneSetCollections, adjMethod = 'fdr',nPerm = 1000, ncpus = 4)
    gsaRes5 <- runGSA(pstat, geneSetStat="wilcoxon", directions = fold_change, gsc=geneSetCollections, adjMethod = 'fdr',nPerm = 1000, ncpus = 4)
    
    #plotting the consensus results
    resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5)
    names(resList) <- c("fisher","stouffer","reporter","tailStrength","wilcoxon")
    ch <- consensusHeatmap(resList,cutoff=50,method="median", ncharLabel = 50, cellnote = "medianPvalue", cex = 0.2, plot = F) 
    
    
    #INTERSECTING TERMS 
    png(filename = "PLOT/NonDirectional.png",width = 9, height = 8)
    cs_nonDirection <- consensusScores(resList, class="non",plot = T)
    dev.off()
    png(filename ="PLOT/DistinctUP.png",width = 9, height = 8)
    cs_Distinct_up <- consensusScores(resList, class="distinct", direction = "up")
    dev.off()
    png(filename ="PLOT/MixedUP.png",width = 9, height = 8)
    cs_Mixed_up <- consensusScores(resList, class="mixed", direction = "up")
    dev.off()
    png(filename ="PLOT/DistinctDOWN.png",width = 9, height = 8)
    cs_Distinct_down <- consensusScores(resList, class="distinct", direction = "down")
    dev.off()
    png(filename ="PLOT/MixedDOWN.png",width = 9, height = 8)
    cs_Mixed_down <- consensusScores(resList, class="mixed", direction = "down")
    dev.off()
    
    if (plots == TRUE){
        
        
        
        # INTERSECT NON DIRECTIONAL AND DIRECTIONAL UP
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_nonDirection$rankMat)), 
                           length(rownames(cs_Distinct_up$rankMat)), 
                           length(intersect(rownames(cs_nonDirection$rankMat), 
                                            rownames(cs_Distinct_up$rankMat))), category = c("NON DIRECTIONAL", "DISTINCT UP"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        
        # INTERSECT NON DIRECTIONAL AND DIRECTIONAL DOWN
        
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_nonDirection$rankMat)), 
                           length(rownames(cs_Distinct_down$rankMat)), 
                           length(intersect(rownames(cs_nonDirection$rankMat), 
                                            rownames(cs_Distinct_down$rankMat))), category = c("NON DIRECTIONAL", "DISTINCT DOWN"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        # INTERSECT NON DIRECTIONAL AND MIXED UP
        
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_nonDirection$rankMat)), 
                           length(rownames(cs_Mixed_up$rankMat)), 
                           length(intersect(rownames(cs_nonDirection$rankMat), 
                                            rownames(cs_Mixed_up$rankMat))), category = c("NON DIRECTIONAL", "MIXED UP"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        # INTERSECT NON DIRECTIONAL AND MIXED DOWN
        
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_nonDirection$rankMat)), 
                           length(rownames(cs_Mixed_down$rankMat)), 
                           length(intersect(rownames(cs_nonDirection$rankMat), 
                                            rownames(cs_Mixed_down$rankMat))), category = c("NON DIRECTIONAL", "MIXED DOWN"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        # INTERSECT DIRECTIONAL UP AND MIXED UP
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_Distinct_up$rankMat)), 
                           length(rownames(cs_Mixed_up$rankMat)), 
                           length(intersect(rownames(cs_Distinct_up$rankMat), 
                                            rownames(cs_Mixed_up$rankMat))), category = c("DIRECTIONAL UP", "MIXED UP"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        
        # INTERSECT DIRECTIONAL DOWN AND MIXED DOWN
        grid.newpage()
        draw.pairwise.venn(length(rownames(cs_Mixed_down$rankMat)), 
                           length(rownames(cs_Distinct_down$rankMat)), 
                           length(intersect(rownames(cs_Mixed_down$rankMat), 
                                            rownames(cs_Distinct_down$rankMat))), category = c("DIRECTIONAL DOWN", "MIXED DOWN"), lty = rep("blank", 2), 
                           fill = c("light blue", "pink"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = rep(0.025, 2))
        
        
        
        
        # INTERSECT AMONG NON DIRECTION, DIRECTIONAL UP AND MIXED UP
        
        png(filename="PLOT/UPIntersection.png")
        draw.triple.venn(area1 = length(rownames(cs_Distinct_up$rankMat)), 
                         area2 = length(rownames(cs_Mixed_up$rankMat)), 
                         area3 = length(rownames(cs_nonDirection$rankMat)), 
                         n12 = length(intersect(rownames(cs_Distinct_up$rankMat), 
                                                rownames(cs_Mixed_up$rankMat))), 
                         n23 = length(intersect(rownames(cs_nonDirection$rankMat), 
                                                rownames(cs_Mixed_up$rankMat))), 
                         n13 = length(intersect(rownames(cs_nonDirection$rankMat), 
                                                rownames(cs_Distinct_up$rankMat))), 
                         n123 = length(intersect(intersect(rownames(cs_Distinct_up$rankMat), 
                                                           rownames(cs_Mixed_up$rankMat)),rownames(cs_nonDirection$rankMat))), 
                         category = c("DIRECTIONAL UP", "MIXED UP", "N0N DIRECTIONAL"), 
                         lty = "blank", 
                         fill = c("skyblue", "pink1", "mediumorchid"))
        dev.off()
       
        # INTERSECT AMONG NON DIRECTION, DIRECTIONAL DOWN AND MIXED DOWN
        
        png(filename="PLOT/DOWNIntersection.png")
        draw.triple.venn(area1 = length(rownames(cs_Distinct_down$rankMat)), 
                         area2 = length(rownames(cs_Mixed_down$rankMat)), 
                         area3 = length(rownames(cs_nonDirection$rankMat)), 
                         n12 = length(intersect(rownames(cs_Distinct_down$rankMat), 
                                                rownames(cs_Mixed_down$rankMat))), 
                         n23 = length(intersect(rownames(cs_nonDirection$rankMat), 
                                                rownames(cs_Mixed_down$rankMat))), 
                         n13 = length(intersect(rownames(cs_nonDirection$rankMat), 
                                                rownames(cs_Distinct_down$rankMat))), 
                         n123 = length(intersect(intersect(rownames(cs_Distinct_down$rankMat), 
                                                           rownames(cs_Mixed_down$rankMat)),rownames(cs_nonDirection$rankMat))), 
                         category = c("DIRECTIONAL DOWN", "MIXED DOWN", "N0N DIRECTIONAL"), 
                         lty = "blank", 
                         fill = c("skyblue", "pink1", "mediumorchid"))
        dev.off()
        
    }
    
   up_bp <- intersect(intersect(rownames(cs_Distinct_up$rankMat), 
                                rownames(cs_Mixed_up$rankMat)),rownames(cs_nonDirection$rankMat))
   
   down_bp <- intersect(intersect(rownames(cs_Distinct_down$rankMat), 
                                  rownames(cs_Mixed_down$rankMat)),rownames(cs_nonDirection$rankMat))
   
   newList <- list("up" = up_bp, "down" = down_bp)
   return (newList)
    
    
}
