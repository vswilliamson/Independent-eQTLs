"plotindividualr2" <-
function(n = 5, li.unipval, li.bipval, li.uniaic, li.biaic)
{
    source("calcanova.R")
    df.unipval <- calcanova(df = li.unipval[[1]], snps = li.unipval[[2]], n = n)
    df.bipval <- calcanova(df = li.bipval[[1]], snps = li.bipval[[2]], n = n)
    df.uniaic <- calcanova(df = li.uniAIC[[1]], snps = li.uniAIC[[2]], n = n)
    df.biaic <- calcanova(df = li.biAIC[[1]], snps = li.biAIC[[2]], n = n)
    
    names <- c("unipval", "bipval", "uniaic", "biaic")
    li <- list(df.unipval, df.bipval, df.uniaic, df.biaic)
    factors <- factor()
    x <- numeric()
    y <- numeric()
    se <- numeric()
    for(i in 1:4){
        m <- length(li[[i]][,1])
        for(j in 1:n){
            factors <- append(factors, names[i])
            x <- append(x, j)
            y <- append(y, mean(li[[i]][,j], na.rm = TRUE))
            se <- append(se, sd(li[[i]][,j])/sqrt(length(li[[i]][,j])))
        }
    }
    return(data.frame(factors, x, y, se))
}

