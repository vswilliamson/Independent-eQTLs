"setUpEnvironment" <-
function(geneExpressionFile, genotypeFile, genelocFile, snpslocFile, me)
{
    GE <- read.table(geneExpressionFile, header = TRUE)
    row.names(GE) <- GE[,1]
    GE <- GE[,-1]
    GE.train <- GE[,1:floor(ncol(GE))]
    GE.test <- GE[,floor(ncol(GE)) + 1:ncol(GE)]

    SNP <- read.table(genotypeFile, header = TRUE, na.strings = c("-1"))
    row.names(SNP) <- SNP[,1]
    row.names(SNP) <- gsub(":", ".", row.names(SNP.t))
    SNP <- SNP[,-1]
    SNP.t <- as.data.frame(t(SNP))
    SNP.train <- SNP[1:floor(nrow(SNP.t)),]
    SNP.test <- SNP[floor(nrow(SNP.t)) + 1:nrow(GE),]

    eQTLs <- me$cis$eqtls
    eQTLs$snps <- gsub(":", ".", eQTLs$snps)
    besteQTLs <- eQTLs[!duplicated(eQTLs$gene),]
    row.names(besteQTLs) <- besteQTLs$gene

    return(list(GE = GE, GE.train = GE.train, GE.test = GE.test, SNP.t = SNP.t, SNP.train = SNP.train, SNP.test = SNP.test, eQTLs = eQTLs, besteQTLs = besteQTLs))
}
