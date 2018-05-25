#Function to perform Quantile normalization
quantileNormalization <-function(wd, distribution) {
    
    n <- nrow(wd)
    m <- ncol(wd)
    if(!missing(distribution)) 
        if(m != length(distribution))
            stop("The reference distribution has length different from the col dimension of the data matrix.") else
                distribution  <-  sort(distribution)
    
    o <- matrix(0, n, m)
    for(i in 1:n)
        o[i,] <- order(wd[i,])
    
    j <- 1
    tmp <- rep(0, n)
    
    while(j <= m) {
        if(missing(distribution)) {
            for(i in 1:n)
                tmp[i] <- wd[i,o[i,j]]
            value <- mean(tmp)
        } else value  <- distribution[j]
        for(i in 1:n)
            wd[i,o[i,j]] <- value
        j <- j + 1
    }
    return(wd)
}