
rm(list=ls()) #clear screen

##set working dir
setwd("C:/Users/dell/Desktop/Data")

##used packages
library(MASS)
library(regDIF)

bt<-proc.time()

##true values
data1 <- read.csv(file="b21f.csv", header = TRUE)
head(data1)
dim(data1)

data11 <- data1[(data1$age>=16)&(data1$age<=25),]
head(data11)
dim(data11)

y <- as.matrix(data11[,c(3:21)])
x <- as.matrix(cbind(data11[,2], sd(data11[,2])*scale(data11[,1])))

N=dim(y)[1]
J=dim(y)[2]
P=dim(x)[2]

##model fit
set.seed(10)

aa <- regDIF(item.data=y, pred.data=x, item.type="2pl", pen.type="lasso", tau=seq(17.5,13.5,length.out=100), stdz=FALSE)

op <- max( which( aa$bic == min(na.omit(aa$bic)) ) )

res <- c(coef(aa)$base[1:J,op], coef(aa)$dif[1:(2*J),op], coef(aa)$base[(J+1):(2*J),op], coef(aa)$dif[(2*J+1):(4*J),op], coef(aa)$impact[,op])
  
et <- proc.time()
 
print((et-bt)[3])

res

save.image(paste("Freq922",".RData",sep=""))











