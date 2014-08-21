"testAll" <-
function(GE, SNPs, eQTLs, method = "pvalue", direction = "forward", steps = 100, verbose = TRUE)
{   
    source("stepwise.R")

    genes <- character()
    numind <- numeric()
    rsquared <- numeric()
    snps <- list()
    i <- 1
    besteQTLs <- eQTLs[!duplicated(eQTLs$gene),]
    for(gene in besteQTLs$gene) {
        gene <- toString(gene)
        if(verbose) cat("Gene #", i, ": ", gene, "\n", sep="")
        tmp <- summary(stepwise(GE = GE, SNPs = SNPs, eQTLs = eQTLs, gene = gene, method = method, direction = direction, steps = steps))
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
        genes <- append(genes, gene)
        numind <- append(numind, num)
        rsquared <- append(rsquared, r2)
        snps[length(snps) + 1] <- list(snp.add)
        i <- i + 1
        if(!verbose) next
        if(is.na(num)){
            cat("0 independent eQTLs; there were no rows in SNPs with all non-NA values.\n")
            next
        }
        cat(num, "independent eQTLs:", paste(snp.add), "\n")
    }
    indeQTLs <- data.frame(genes, numindeQTLs=numind, rsquared)
    if(verbose) cat("Done.\n") 
    return(list(indeQTLs=indeQTLs, snps=snps))
}
