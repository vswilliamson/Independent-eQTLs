"crossvalidate" <-
function(eQTLs, sortBy = c("beta", "pvalue"))
{
    source("stepwise.R")
    methods <- c("pvalue", "AIC")
    directions <- c("forward", "both")
    firstmean <- matrix(ncol=2, nrow=2)
    secondmean <- matrix(ncol=2, nrow=2)
    trainr2 <- list(list(numeric(), numeric()), list(numeric(), numeric()))
    testr2 <- list(list(numeric(), numeric()), list(numeric(), numeric()))
    source("r2.R")

    numeqtl <- numeric()
    genes <- numeric()
    r2 <- numeric()
    method <- character()
    set <- character()
    mp <- numeric()
    for(gene in besteQTLs$gene) {
        gene <- toString(gene)
        has.na <- FALSE
        train <- list(list(), list())
        for(i in 1:2){
            for(j in 1:2){  
                cat(methods[i], directions[j], gene, "\n")
                tmp <- stepwise(GE = GE.train, SNPs = SNP.train, eQTLs, gene = gene, method = methods[i], direction = directions[j])
                train[[i]][[j]] <- tmp
                addtrain <- NULL
                addtest <- NULL
                if(is.null(summary(tmp)) || summary(tmp)[2] == "NULL") {
                    addtrain <- NA
                    addtest <- NA
                    has.na <- TRUE
                }
                else {
                    addtrain <- summary(tmp)$r.squared
                    addtest <- r2(f=as.numeric(predict.lm(object=tmp, newdata=SNP.test)), y=as.numeric(GE.test[gene,]))
                }
                trainr2[[i]][[j]] <- append(trainr2[[i]][[j]], addtrain)
                testr2[[i]][[j]] <- append(testr2[[i]][[j]], addtest)
            }
        }
        m <- length(trainr2[[1]][[1]])
        for(i in 1:2){
            for(j in 1:2){
                trainr2[[i]][[j]] <- trainr2[[i]][[j]][-m]
                testr2[[i]][[j]] <- testr2[[i]][[j]][-m]
            }
        }
        if(!has.na) {
            minpredictors <- 1000000 #arbitrarily big number
            for(i in 1:2){
                for(j in 1:2){
                    minpredictors <- min(minpredictors, train[[i]][[j]]$rank - 1)
                }
            }
            mp <- append(mp, minpredictors)
            for(i in 1:2){
                for(j in 1:2){
                    lm.train <- train[[i]][[j]]
                    
                    lm.train <- summary(lm.train)$coefficients
                    if(nrow(lm.train) == 2) {
                        second.name <- row.names(lm.train)[2]
                    }
                    lm.train <- lm.train[-1,]
                    if(is.null(dim(lm.train))) {
                        lm.train <- as.data.frame(rbind(lm.train))
                        row.names(lm.train) <- second.name
                    }
                    names.train <- NULL
                    if(sortBy == "beta") {
                        names.train <- row.names(lm.train[order(-abs(lm.train[,1])),])[1:minpredictors]
                    }
                    else {
                        names.train <- row.names(lm.train[order(lm.train[,4]),])[1:minpredictors]
                    }
                    y.tmp <- as.numeric(GE.train[gene,])
                    lm.train <- lm(formula = y.tmp~.,data = SNP.train[names.train])
                    addtrain <- summary(lm.train)$r.squared
                    tmp <- na.omit(data.frame(as.numeric(predict.lm(object=lm.train, newdata=SNP.test)), as.numeric(GE.test[gene,])))
                    addtest <- (cor(tmp[,1], tmp[,2]))^2
                    trainr2[[i]][[j]] <- append(trainr2[[i]][[j]], addtrain)
                    testr2[[i]][[j]] <- append(testr2[[i]][[j]], addtest)
                    genes <- append(genes, gene)
                    genes <- append(genes, gene)
                    method <- append(method, paste(methods[i], directions[j]))
                    method <- append(method, paste(methods[i], directions[j]))
                    set <- append(set, "train")
                    set <- append(set, "test")
                    numeqtl <- append(numeqtl, minpredictors)
                    numeqtl <- append(numeqtl, minpredictors)
                    r2 <- append(r2, addtrain)
                    r2 <- append(r2, addtest)
                }
            }
        }
    }
    return(data.frame(genes, method, set, numeqtl, r2))
}
