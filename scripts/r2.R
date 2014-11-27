"r2" <-
function (y, f)
{
    mat <- cbind(y, f)
    mat <- na.omit(mat)
    y <- as.numeric(mat[,1])
    f <- as.numeric(mat[,2])
    1-(sum((y-f)^2, na.rm=TRUE)/sum((y-mean(y))^2, na.rm=TRUE))
}
