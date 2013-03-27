####checks for package presence and load it
if (length(grep('grofit',installed.packages()[,1]))==0)
{
 install.packages("grofit")
}
library(grofit)

####function

analyse_growth=function(file,export=T,name,y_scale=c(0.8,3.3))
{
 
 library(grofit)
 data=read.delim(file,stringsAsFactors=F)
 out=data.frame(stringsAsFactors=F)
 if(export==T)
 {
 old=getwd()
 dir.create(paste("growth_analysis_",name,sep=''))
 setwd(paste("growth_analysis_",name,sep=''))
 }
 for(i in 2:ncol(data))
 #for(i in 3:4)
 {
###
print(i-1)
#i=18
  j=i-1
  fit1=gcFitSpline(time=data[,1], data=data[,i])
  out[j,1]=colnames(data)[i]
  out[j,2]=fit1$parameters$A
  out[j,3]=fit1$parameters$lambda
  out[j,4]=fit1$parameters$mu
  out[j,5]=fit1$parameters$integral
#print(i)
  inf1=fit1$fit.data[2:length(fit1$fit.data)]
  inf2=fit1$fit.data[1:length(fit1$fit.data)-1]
  inf3=inf1-inf2
  inf4=cbind(fit1$fit.time,c(0,inf3),fit1$fit.data)
  inf5=inf4[order(inf4[,2]),]
  inf6=inf5
#print(inf6)
  inf5=inf5[which(inf5[,1] > inf5[nrow(inf5),1]),]
##
#print(i)
#print(inf5)
#print(nrow(inf5))
##
  #if(nrow(inf5) > 0)
  if((is.null(nrow(inf5))==FALSE)) 
  {
   if(nrow(inf5) > 0)
   {
    out[j,6]=inf5[1,3]
    out[j,7]=inf5[2,3]
    out[j,8]=inf5[3,3]
    out[j,9]=inf6[nrow(inf6),1]
   }
   else
   {
    out[j,6]=NA
    out[j,7]=NA
    out[j,8]=NA
    out[j,9]=inf6[nrow(inf6),1]
   }
  }
  else
  {
   out[j,6]=NA
   out[j,7]=NA
   out[j,8]=NA
   out[j,9]=inf6[nrow(inf6),1]
  }


#print(i)
###
#setwd("C:/Documents and Settings/Samuel/My Documents/")
#return(fit1)
###
 
  if(export==T)
  {
   jpeg(file=paste("picture_",name,"_",j,".jpg",sep=''))
   plot.growth(fit1,main=out[j,1],scale=y_scale)
   abline(h=fit1$parameters$A,col="red",lwd=2)
   abline(v=fit1$parameters$lambda,col="red",lwd=2)
   abline(h=out[j,7],col="grey",lwd=2)
   abline(h=out[j,8],col="grey",lwd=2)
   abline(h=out[j,6],col="blue",lwd=2)
   abline(v=out[j,9],col="green",lwd=2)

   legend(x="bottomright",legend=colnames(data)[i],cex=1.5,bty="n",box.col="white",bg="white")
   dev.off()
  }
  rm(fit1)
 }
 colnames(out)=c("name","max.mass","lag","max.slope","integral","Plateau.1","Plateau.2","Plateau.3","max.slope.time")

 if (export==T)
 {
  write.table(out,file=paste("growth_data_",name,".txt",sep=''),row.names=F,quote=F,sep="\t")
  setwd(old)
 }

return(data)

}
##################################################################
plot.growth=function (x, add = FALSE, raw = TRUE, slope = TRUE, pch = 1, 
    colData = 1, colSpline = 2, cex = 1, scale=c(0.8,3.3),...) 
{
    if (is.logical(add) == FALSE) 
        stop("Need logical value for: add")
    if (is.logical(raw) == FALSE) 
        stop("Need logical value for: raw")
    if (is.logical(slope) == FALSE) 
        stop("Need logical value for: slope")
    if (is.numeric(pch) == FALSE) 
        stop("Need numeric value for: pch")
    if (is.numeric(cex) == FALSE) 
        stop("Need numeric value for: cex")
    if (FALSE %in% (colData %in% c(colors(), 0:8))) 
        stop("col needs to be numeric from 0:8 or a string from colors()")
    if (FALSE %in% (colSpline %in% c(colors(), 0:8))) 
        stop("col needs to be numeric from 0:8 or a string from colors()")
    if ((is.na(x$fitFlag) == TRUE) | (x$fitFlag == FALSE)) {
        warning("plot.gcFitModel: no data fit available!")
    }
    else {
        if (raw == TRUE) {
            if (add == TRUE) {
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(points(x$raw.time, x$raw.data, sub = x$name.fit, 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(points(x$raw.time, x$raw.data, sub = x$name.fit, 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(points(x$raw.time, x$raw.data, sub = x$name.fit, 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(points(x$raw.time, x$raw.data, sub = x$name.fit, 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
            }
            else {
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  FALSE)) {
                  if(is.null(scale)==TRUE)
                  {
                  try(plot(x$raw.time, x$raw.data, sub = x$name.fit, 
                    xlab = "time", ylab = "growth y(t)", col = colData, 
                    pch = pch, cex = cex))
                  }
                  else
                  {
                  try(plot(x$raw.time, x$raw.data, sub = x$name.fit, 
                    xlab = "time", ylab = "growth y(t)", col = colData, 
                    pch = pch, cex = cex, ylim=c(scale[1],scale[2])))
                  }
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(plot(x$raw.time, x$raw.data, sub = x$name.fit, 
                    xlab = "time", ylab = "log(1+growth y(t))", 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(plot(x$raw.time, x$raw.data, sub = x$name.fit, 
                    xlab = "log(1+time)", ylab = "growth y(t)", 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(plot(x$raw.time, x$raw.data, sub = x$name.fit, 
                    xlab = "log(1+time)", ylab = "log(1+growth y(t))", 
                    col = colData, pch = pch, cex = cex))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
            }
        }
        else {
            if (add == TRUE) {
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
            }
            else {
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(plot(x$fit.time, x$fit.data, sub = x$name.fit, 
                    xlab = "time", ylab = "growth y(t)", type = "n"))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == FALSE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(plot(x$fit.time, x$fit.data, sub = x$name.fit, 
                    xlab = "time", ylab = "log(1+growth y(t))", 
                    type = "n"))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  FALSE)) {
                  try(plot(x$fit.time, x$fit.data, sub = x$name.fit, 
                    xlab = "log(1+time)", ylab = "growth y(t)", 
                    type = "n"))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
                if ((x$control$log.x.gc == TRUE) && (x$control$log.y.gc == 
                  TRUE)) {
                  try(plot(x$fit.time, x$fit.data, sub = x$name.fit, 
                    xlab = "log(1+time)", ylab = "log(1+growth y(t))", 
                    type = "n"))
                  try(lines(x$fit.time, x$fit.data, sub = x$name.fit, 
                    col = colSpline, type = "l"))
                }
            }
        }
        if (slope == TRUE) {
            mu <- as.numeric(x$parameters$mu)
            lambda <- as.numeric(x$parameters$lambda)
            bla <- x$fit.time * mu
            bla <- bla + (-x$parameters$mu * x$parameters$lambda)
            try(lines(x$fit.time, bla, lw = 2, lty = 2, col = colSpline))
        }
    }
}

