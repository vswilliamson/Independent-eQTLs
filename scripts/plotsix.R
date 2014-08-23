"plotsix" <-
#6x1 plot, mar: margin, oma: outer margin
function(df, snps, sortBy="beta")
{
    if(sortBy != "beta" && sortBy != "pvalue") stop("sortBy must equal \"beta\" or \"pvalue\".")
    source("firstfive.R")
    li.tmp <- firstfive(df, snps, sortBy=sortBy)
    par(mfrow=c(6,1), mar=c(1,3,1,1), oma=c(1.5,0,0.1,0))

    plot(li.tmp[[12]], abs(li.tmp[[11]]), ylim=c(0, 1), xlim=c(-1000000, 1000000), main="Best eQTL reported by Matrix eQTL", xlab="", ylab="", sub = "", xaxt='n')

    for(i in 1:4)
    {
        plot(li.tmp[[i*2]], abs(li.tmp[[i*2-1]]), ylim=c(0, 1), xlim=c(-1000000, 1000000), main=paste("SNP #", i, sep=""), xlab="", ylab="", sub = "", xaxt='n')
    }   

    plot(li.tmp[[10]], abs(li.tmp[[9]]), ylim=c(0, 1), xlim=c(-1000000, 1000000), main="SNP #5", xlab="", ylab="Effect size", sub = "")
    print(sd(li.tmp[[12]]))
    for(i in 1:5)
    {
        print(sd(li.tmp[[i*2]]))
    }
}
