#Adapted from http://www.bioconductor.org/packages/release/bioc/html/maSigPro.html

"stepwise" <-
function (GE, SNPs, eQTLs, gene, method = "pvalue", direction = "forward", steps = 100)
{
    y <- as.numeric(GE[gene,])
    d <- SNPs[as.character(eQTLs[eQTLs$gene == gene, 1])]
    if(nrow(na.omit(d)) == 0)
    { 
        return(NULL)
    }

    ##AIC selection
    if (method == "AIC") {
        if(nrow(summary(lm(y ~ ., data = d, na.action = na.exclude))$coefficients) == 1)
        {
            return(NULL)
        }
        resul0 <- summary(lm(y ~ ., data = d, na.action = na.exclude))$coefficients[, 4]
        tmp <- gsub("`", "", names(resul0)[-1])
        d <- as.data.frame(d[tmp])
        tmp <- cbind(d, y)
        tmp <- na.omit(tmp)
        y.no.na <- tmp$y
        d.no.na <- tmp
        d.no.na$y <- NULL
        null <- lm(formula = y.no.na~1, data = d.no.na)
        full <- lm(formula = y.no.na~., data = d.no.na)
        return(step(null, scope = list(lower = null, upper = full), direction = direction, steps = steps, trace = 0))
    }

    #DEPRECIATED: lars lasso
    if (method == "lars") {
        library("lars")
        if(nrow(summary(lm(y ~ ., data = d, na.action = na.exclude))$coefficients) == 1)
        {
          return(NULL)
        }
        resul0 <- summary(lm(y ~ ., data = d, na.action = na.exclude))$coefficients[, 4]
        tmp <- gsub("`", "", names(resul0)[-1])
        d <- as.data.frame(d[tmp])
        tmp <- cbind(d, y)
        tmp <- na.omit(tmp)
        y.no.na <- tmp$y
        d.no.na <- tmp
        d.no.na$y <- NULL
        return(lars(x = as.matrix(d.no.na), y = y.no.na, normalize=TRUE, type="lasso"))
    }
    
    #pvalue selection
    alphain <- 0.05;
    alphaout <- 0.05;
    pval <- NULL;
    design <- NULL;
    j = 1;
    resul0 <- summary(lm(y ~ ., data = d, na.action = na.exclude))$coefficients[, 4] #gets rid of collinear variables
    tmp <- gsub("`", "", names(resul0)[-1])
    if(length(tmp) != 0)
    {
        d <- as.data.frame(d[tmp])
    }
    if(ncol(d) == 1) {
        test <- lm(y ~ ., data = d, na.action = na.exclude)
        if( summary(test)$coefficients[,4][2] < 0.05) {
            return(test)  
        }
        else {
            return(NULL)
        }
    }
    for (i in 1:ncol(d)) {
        sub <- cbind(design, d[, i])
        sub <- as.data.frame(sub)
        lm2 <- lm(y ~ ., data = sub, na.action = na.exclude)
        result <- summary(lm2)
        pval[i] <- result$coefficients[, 4][j + 1]
    }
    min <- min(pval, na.rm = TRUE)
    while (min <= alphain && steps > 0 ) {
        steps <- steps - 1
        b <- pval == min
        c <- c(1:length(pval))
        pos <- c[b]
        pos <- pos[!is.na(pos)][1]
        design <- cbind(design, d[, pos])
        design <- as.data.frame(design)
        colnames(design)[j] <- colnames(d)[pos]
        if (ncol(d) == 2) {
            lastname <- colnames(d)[colnames(d) != colnames(d)[pos]]
        }
        d <- d[, -pos]
        if (is.null(dim(d))) {
            d <- as.data.frame(d)
            colnames(d) <- lastname
        }
        result2 <- summary(lm(y ~ ., data = design, na.action = na.exclude))$coefficients[,4]
        max <- max(result2[-1], na.rm = TRUE)
        if(direction != "forward")
        {
            while (max > alphaout) {
                varout <- names(result2)[result2 == max]
                varout <- gsub("`", "", varout)
                pos <- position(matrix = design, vari = varout)
                d <- as.data.frame(cbind(d, design[, pos]))
                x <- ncol(d)
                colnames(d)[x] <- colnames(design)[pos]
                if (ncol(design) == 2) {
                    min <- min(result2[-1], na.rm = TRUE)
                    lastname <- names(result2)[result2 == min]
                }
                design <- design[, -pos]
                if (is.null(dim(design))) {
                    design <- as.data.frame(design)
                    colnames(design) <- lastname
                 }
                result2 <- summary(lm(y ~ ., data = design, na.action = na.exclude))$coefficients[,4]
                max <- max(result2[-1], na.rm = TRUE)
            }
        }
        j = ncol(design) + 1
        pval <- NULL
        for (i in 1:ncol(d)) {
            sub <- cbind(design, d[, i])
            sub <- as.data.frame(sub)
            lm2 <- lm(y ~ ., data = sub , na.action = na.exclude)
            result <- summary(lm2)
            pval[i] <- result$coefficients[, 4][j + 1]
        }
        min <- min(pval, na.rm = TRUE)
        if (ncol(d) == 1) {
            if (min <= alphain) {
                design <- cbind(design, d[, 1])
                design <- as.data.frame(design)
                colnames(design)[j] <- colnames(d)[1]
            }
            min = 1
        }
    }
    if (is.null(design)) {
        lm1 <- lm(y ~ 1, na.action = na.exclude)
    }
    else {
        lm1 <- lm(y ~ ., data = design, na.action = na.exclude)
    }
    names(lm1$coefficients) <- gsub("`", "", gsub("\\\\", "", names(lm1$coefficients)))
    return(lm1)
}

