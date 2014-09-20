"calcanova" <-
function(df, snps, n = 5)
{
    ret <- data.frame(matrix(ncol=n, nrow=0))
    for(i in 1:nrow(histo))
    {
        if(!is.na(histo[i, 2])){
            gene <- toString(histo[i,1])
            #if(gene != "ACOT8") next
            y.tmp <- as.numeric(GE[gene,])
            d <- as.data.frame(SNP.t[as.character(snps[[i]])])
            #if(ncol(d) == 1){
            #    colnames(d)[1] <- as.character(snps[[i]])
            #}
            fit <- lm(formula = y.tmp~., data = d)
            af <- anova(fit)
            afss <- af$"Sum Sq"
            propExp <- afss/sum(afss)
            to.sort <- summary(fit)$coefficients[-1,]
            if(ncol(as.data.frame(to.sort)) == 1) next
            to.sort <- cbind(to.sort, propExp=propExp[-length(propExp)])
            to.sort <- to.sort[order(abs(to.sort[,4])),]
            if(nrow(to.sort) >= n) {
                row.names(to.sort) <- gsub("`", "", row.names(to.sort))
                to.add <- unname(to.sort[1:n,5])
                ret <- rbind(ret, to.add)
            }
        }
    }
    return(ret)
    #for(i in 1:5) print(mean(df[,i]))
}

