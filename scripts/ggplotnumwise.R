"ggplotnumwise" <-
function(df, directory, n = 5, stratify = TRUE)
{
    library("ggplot2")
    source("multiplot.R")
    methods = c("pvalue forward", "pvalue both", "AIC forward", "AIC both")
    for(i in 1:n) {
        df.tmp <- df[df$numeqtl == i,]
        if(!stratify) {
            png(paste(directory, i, ".png", sep = ""), 900, 900)
            p <- ggplot(df.tmp, aes(x=r2, fill=set)) + geom_density(alpha=.3) + xlim(0,1) + ggtitle(paste(i, "independent eQTLs"))
            print(p)
            dev.off()
        }
        smalldfs <- lapply(1:4, function(x){df.tmp[df.tmp$method == methods[x],]})
        p <- ggplot(smalldfs[[1]], aes(x=r2, colour=set)) + geom_density(alpha=.3) + xlim(0,1)
        plots <- lapply(smalldfs, FUN = function(x,p){p %+% x},p=p)
        for(j in 1:4){
            plots[[j]] <- plots[[j]] + ggtitle(methods[j])
        }
        setEPS()
        postscript(paste(directory, i, ".eps", sep = ""))
        multiplot(plotlist = plots, cols = 2, main = paste(i, "independent eQTLs"))
        dev.off()
    }
}
