####checks for package presence and load it
if (length(grep('grofit',installed.packages()[,1]))==0)
{
 install.packages("grofit")
}
library(grofit)

####function

analyse_growth=function(file,export=c("short","all","none"),name,y_scale=c(0.8,3.3),discard_up=NULL)
{
 
 library(grofit)
 data=read.delim(file,stringsAsFactors=F)
 out=data.frame(stringsAsFactors=F)
 export=export[1]
 if(export != "none")
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
  if(is.null(discard_up)==FALSE)
  {
   data[which(data[,1] < discard_up),i]=1
  }
  fit1=gcFitSpline(time=data[,1], data=data[,i])
  out[j,1]=colnames(data)[i]
  out[j,2]=fit1$parameters$A
  out[j,3]=fit1$parameters$lambda
  out[j,4]=fit1$parameters$mu
  out[j,5]=fit1$parameters$integral


######extract mass at first fitted point after lag time
#print("ok.1")
TLag=min(fit1$fit.time[(fit1$fit.time-out[j,3])>0])
MassLag=fit1$fit.data[(fit1$fit.time==TLag)]
#print(out[j,3])
#print(TLag)
#print(MassLag)
######


#print("ok.2")
#print(fit1$fit.time)
##
#print(i)
  inf1=fit1$fit.data[2:length(fit1$fit.data)]
  inf2=fit1$fit.data[1:length(fit1$fit.data)-1]
  inf3=inf1-inf2
  inf4=cbind(fit1$fit.time,c(0,inf3),fit1$fit.data)
  inf5=inf4[order(inf4[,2]),]
  inf6=inf5
  inf5=inf5[which(inf5[,1] > inf5[nrow(inf5),1]),]
#print(inf5)
  if(nrow(inf5) >= 20)
  {
   inf7=inf5[1:20,]
   inf7=inf7[which(inf7[,1]==min(inf7[,1])),3]
  }
  else if (nrow(inf5) != 0)
  {
   inf7=inf5[1:nrow(inf5),]
   inf7=inf7[which(inf7[,1]==min(inf7[,1])),3]
  }
  else
  {
   inf7=NA 
  }

#print("ok.3")
##
#print(i)
#print(inf5)
#print(nrow(inf5))
##
  #if(nrow(inf5) > 0)
  if((is.null(nrow(inf5))==FALSE)) 
  {
   if(nrow(inf5) >= 3)
   {
    out[j,6]=inf5[1,3]
    out[j,7]=inf5[2,3]
    out[j,8]=inf5[3,3]
    out[j,9]=inf7
    out[j,10]=inf6[nrow(inf6),1]

TP0=fit1$fit.time[which(fit1$fit.data==out[j,6])]
TP1=fit1$fit.time[which(fit1$fit.data==out[j,7])]
TP2=fit1$fit.time[which(fit1$fit.data==out[j,8])]
TP3=fit1$fit.time[which(fit1$fit.data==out[j,9])]
out[j,11]=(out[j,6]-MassLag)/(TP0-out[j,3])
out[j,12]=(out[j,7]-MassLag)/(TP1-out[j,3])
out[j,13]=(out[j,8]-MassLag)/(TP2-out[j,3])
out[j,14]=(out[j,9]-MassLag)/(TP3-out[j,3])
out[j,15]=TP3
   }
   else
   {
    out[j,6]=NA
    out[j,7]=NA
    out[j,8]=NA
    out[j,9]=NA
    out[j,10]=inf6[nrow(inf6),1]
    out[j,11]=NA
    out[j,12]=NA
    out[j,13]=NA
    out[j,14]=NA
    out[j,15]=NA
   }
  }
  else
  {
   out[j,6]=NA
   out[j,7]=NA
   out[j,8]=NA
   out[j,9]=NA
   out[j,10]=inf6[nrow(inf6),1]
   out[j,11]=NA
   out[j,12]=NA
   out[j,13]=NA
   out[j,14]=NA
   out[j,15]=NA
  }

#print("ok.5")
#print(Plateau0)
#print(Plateau1)
#print(Plateau2)
#print(Plateau3)

#print(i)
###
#setwd("C:/Documents and Settings/Samuel/My Documents/")
#return(fit1)
###
 
  if(export != "none")
  {
   jpeg(file=paste("picture_",name,"_",j,".jpg",sep=''))
   plot.growth(fit1,main=out[j,1],scale=y_scale)
   abline(h=fit1$parameters$A,col="red",lwd=2)
   abline(v=fit1$parameters$lambda,col="red",lwd=2)
   abline(h=out[j,7],col="grey",lwd=2)
   abline(h=out[j,8],col="grey",lwd=2)
   abline(h=out[j,6],col="grey",lwd=2)
   abline(v=out[j,10],col="green",lwd=2)
   abline(h=out[j,9],col="black",lwd=2,lty=2)
   abline(v=out[j,15],col="black",lwd=2,lty=2)
##
abline(h=MassLag,col="black",lwd=2,lty=2)
##
   legend(x="bottomright",legend=colnames(data)[i],cex=1.5,bty="n",box.col="white",bg="white")
   dev.off()
  }
  rm(fit1)
 }
 colnames(out)=c("name","max.mass","lag","max.slope","integral","Plateau.1","Plateau.2","Plateau.3","Plateau.best","max.slope.time","dM/dT.1","dM/dT.2","dM/dT.3","dM/dT.best","t.best")

out[2:15]=round(out[,2:15],2)


 if (export != "none")
 {
  if(export == "all")
  {
  write.table(out,file=paste("growth_data_",name,".txt",sep=''),row.names=F,quote=F,sep="\t")
  }
  if(export == "short")
  {
  write.table(out[,c(1,2,3,4,5,9,10,14)],file=paste("growth_data_",name,".txt",sep=''),row.names=F,quote=F,sep="\t")
  }
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

