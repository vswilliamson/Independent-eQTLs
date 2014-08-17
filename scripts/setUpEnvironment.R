"setUpEnvironment" <-
function(geneExpressionFile, genotypeFile, genelocFile, snpslocFile, me)
{
    GE <- read.table(geneExpressionFile, header = TRUE)
    row.names(GE) <- GE[,1]
    GE <- GE[,-1]
    GE.train <- GE[,1:floor(ncol(GE)/2)]
    GE.test <- GE[,(floor(ncol(GE)/2) + 1):ncol(GE)]

    SNP <- read.table(genotypeFile, header = TRUE, na.strings = c("-1"))
    row.names(SNP) <- SNP[,1]
    row.names(SNP) <- gsub(":", ".", row.names(SNP))
    SNP <- SNP[,-1]
    SNP.t <- as.data.frame(t(SNP))
    SNP.train <- SNP.t[1:floor(nrow(SNP.t)/2),]
    SNP.test <- SNP.t[(floor(nrow(SNP.t)/2) + 1):nrow(SNP.t),]

    eQTLs <- me$cis$eqtls
    eQTLs$snps <- gsub(":", ".", eQTLs$snps)
    besteQTLs <- eQTLs[!duplicated(eQTLs$gene),]
    row.names(besteQTLs) <- besteQTLs$gene

    return(list(GE = GE, GE.train = GE.train, GE.test = GE.test, SNP.t = SNP.t, SNP.train = SNP.train, SNP.test = SNP.test, eQTLs = eQTLs, besteQTLs = besteQTLs))
}
