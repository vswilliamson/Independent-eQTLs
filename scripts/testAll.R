"testAll" <-
function(y, d, method = "pvalue", direction = "forward", steps = 100)
{
    genes <- character()
    numind <- numeric()
    rsquared <- numeric()
    snps <- list()
    histo <- NULL
    i <- 1
    source("stepwise.R")

    for(gene in besteQTLs$gene) {
        gene <- toString(gene)
        print(i)
        print(gene)
        tmp <- summary(stepwise(gene = gene, method = method, direction = direction, steps = steps))
        num <- NULL
        snp.add <- NULL
        if(is.null(tmp) || tmp[2] == "NULL") {
            num <- NA
            r2 <- NA
            snp.add <- NA
        }
        else {
            num <- nrow(tmp$coefficients) - 1
            r2 <- tmp$r.squared
            snp.add <- gsub("`", "", row.names((tmp$coefficients))[-1])
        }
        print(num)
        genes <- append(genes, gene)
        numind <- append(numind, num)
        rsquared <- append(rsquared, r2)
        print(snp.add)
        snps[length(snps) + 1] <- list(snp.add)
        i <- i + 1
    }
    histo <- data.frame(genes, numind, rsquared)
    return(list(histo, snps))
}
