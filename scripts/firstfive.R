"firstfive" <-
function(df, snps, n = 5, sortBy = "beta")
{

    li <- list()

    for(i in 1:12)
    {
        li[[i]] <- numeric()
    }
    for(i in 1:nrow(df))
    {
        if(!is.na(df[i, 2])){
            gene <- toString(df[i,1])
            #if(gene != "ACOT8") next
            y.tmp <- as.numeric(GE[gene,])
            d <- as.data.frame(SNP.t[as.character(snps[[i]])])
            #if(ncol(d) == 1){
            #    colnames(d)[1] <- as.character(snps[[i]])
            #}
            to.sort <- summary(lm(formula = y.tmp~., data = d))$coefficients[-1,]
            if(ncol(as.data.frame(to.sort)) == 1) next
            to.sort <- to.sort[order(-abs(to.sort[,1])),]
            if(sortBy == "pvalue"){
                to.sort <- to.sort[order(to.sort[,4]),]
            }
            if(nrow(to.sort) > n - 1) {
                row.names(to.sort) <- gsub("`", "", row.names(to.sort))
                for(i in 0:4)
                {
                    li[[i*2 + 1]] <- c(li[[i*2 + 1]], to.sort[i+1,1])
                    li[[i*2 + 2]] <- c(li[[i*2 + 2]], geneloc[gene,3] - snpsloc[as.character(row.names(to.sort)[i+1]), 3])
                }
                snp <- as.character(ttest[gene, 1])
                d <- as.data.frame(SNP.t[,snp])
                coef <- summary(lm(formula = y.tmp~., data = d))$coefficients[,1][-1]
                diff <- geneloc[gene,3] - snpsloc[snp, 3]
                li[[2*n + 1]] <- c(li[[2*n + 1]], coef)
                li[[2*n + 2]] <- c(li[[2*n + 2]], diff)
            }
        }
    }
    return(li)
}
