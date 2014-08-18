"setUpEnvironment" <-
function(geneExpressionFile, genotypeFile, genelocFile, snpslocFile, me)
{
    GE <- read.table(geneExpressionFile, header = TRUE)
    GE.train <- GE[,1:floor(ncol(GE)/2)]
    GE.test <- GE[,(floor(ncol(GE)/2) + 1):(2 * floor(ncol(GE)/2))]
    cat("Gene expression file parsed.\n")

    SNP <- read.table(genotypeFile, header = TRUE, na.strings = c("-1"))
    SNP.t <- as.data.frame(t(SNP))
    SNP.train <- SNP.t[1:floor(nrow(SNP.t)/2),]
    SNP.test <- SNP.t[(floor(nrow(SNP.t)/2) + 1):(2 * floor(nrow(SNP.t)/2)),]
    cat("Genotype file parsed.\n")

    geneloc <- read.table(genelocFile, header = TRUE)
    geneloc <- geneloc[!duplicated(geneloc[,1]),]
    row.names(geneloc) <- geneloc[,1]
    cat("Gene location file parsed.\n")

    snpsloc <- read.table(snpslocFile, header = TRUE)
    row.names(snpsloc) <- snpsloc[,1]
    cat("SNP location file parsed.\n")
    
    eQTLs <- me$cis$eqtls
    eQTLs$snps <- gsub(":", ".", eQTLs$snps)
    eQTLs <- eQTLs[order(eQTLs$gene, eQTLs$pvalue),]
    besteQTLs <- eQTLs[!duplicated(eQTLs$gene),]
    row.names(besteQTLs) <- besteQTLs$gene
    cat("Matrix eQTL object parsed.\n")

    return(list(GE = GE, GE.train = GE.train, GE.test = GE.test, SNP.t = SNP.t, SNP.train = SNP.train, SNP.test = SNP.test, geneloc = geneloc, snpsloc = snpsloc, eQTLs = eQTLs, besteQTLs = besteQTLs))
}
