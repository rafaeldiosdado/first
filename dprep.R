mmnorm <- function (data,minval=0,maxval=1) 
{
#This is a function to apply min-max normalization to a matrix or dataframe.
#Min-max normalization subtracts the minimum of an attribute from each value
#of the attribute and then divides the difference by the range of the attribute.
#These new values are multiplied by the given range of the attribute
#and finally added to the given minimum value of the attribute.
#These operations transform the data into [minval,mxval].
#Usually minval=0 and maxval=1.
#Uses the scale function found in the R base package.
#Input: data= The matrix or dataframe to be scaled


#store all attributes of the original data
d=dim(data)
c=class(data)
cnames=colnames(data)

#remove classes from dataset
classes=data[,d[2]]
data=data[,-d[2]]

minvect=apply(data,2,min)
maxvect=apply(data,2,max)
rangevect=maxvect-minvect
zdata=scale(data,center=minvect,scale=rangevect)

#remove attributes added by the function scale and turn resulting
#vector back into a matrix with original dimensions
#attributes(zdata)=NULL
#zdata=matrix(zdata,dim(data)[1],dim(data)[2])

newminvect=rep(minval,d[2]-1)
newmaxvect=rep(maxval,d[2]-1)
newrangevect=newmaxvect-newminvect
zdata2=scale(zdata,center=FALSE,scale=(1/newrangevect))
zdata3=zdata2+newminvect

zdata3=cbind(zdata3,classes)

if (c=="data.frame") zdata3=as.data.frame(zdata3)
colnames(zdata3)=cnames
return(zdata3)

}

signorm<-function (data) 
{
  d = dim(data)
  c = class(data)
  cnames = colnames(data)
  classes = data[, d[2]]
  zdata = znorm(data)
  d2 = dim(zdata)
  zdata = zdata[, -d2[2]]
  sigdata = (1 - exp(-zdata))/(1 + exp(-zdata))
  sigdata = cbind(sigdata, classes)
  if (c == "data.frame") 
    sigdata = as.data.frame(sigdata)
  colnames(sigdata) = cnames
  return(sigdata)
}

decscale<- function (data) 
{
  d = dim(data)
  c = class(data)
  cnames = colnames(data)
  classes = data[, d[2]]
  data = data[, -d[2]]
  maxvect = apply(abs(data), 2, max)
  kvector = ceiling(log10(maxvect))
  scalefactor = 10^kvector
  decdata = scale(data, center = FALSE, scale = scalefactor)
  attributes(decdata) = NULL
  decdata = matrix(decdata, dim(data)[1], dim(data)[2])
  decdata = cbind(decdata, classes)
  if (c == "data.frame") 
    decdata = as.data.frame(decdata)
  colnames(decdata) = cnames
  return(decdata)
}

signorm <- function (data) 
{
  d = dim(data)
  c = class(data)
  cnames = colnames(data)
  classes = data[, d[2]]
  zdata = znorm(data)
  d2 = dim(zdata)
  zdata = zdata[, -d2[2]]
  sigdata = (1 - exp(-zdata))/(1 + exp(-zdata))
  sigdata = cbind(sigdata, classes)
  if (c == "data.frame") 
    sigdata = as.data.frame(sigdata)
  colnames(sigdata) = cnames
  return(sigdata)
}

znorm <-function (data) 
{
  d = dim(data)
  c = class(data)
  cnames = colnames(data)
  classes = data[, d[2]]
  data = data[, -d[2]]
  zdata = scale(data)
  attributes(zdata) = NULL
  zdata = matrix(zdata, dim(data)[1], dim(data)[2])
  zdata = cbind(zdata, classes)
  if (c == "data.frame") 
    zdata = as.data.frame(zdata)
  colnames(zdata) = cnames
  return(zdata)
}

disc.ef<-function (data, varcon, k) 
{
  p <- dim(data)[2]
  f <- p - 1
  ft <- rep(0, f)
  for (i in 1:length(varcon)) {
    ft[varcon[i]] = 1
  }
  for (i in 1:f) {
    if (ft[i] > 0) {
      data[, i] <- disc2(as.vector(data[, i]), k)
    }
  }
  data
}

EqualFreq2 <- function(x,n){
  nx <- length(x)
  nrepl <- floor(nx/n)
  nplus <- sample(1:n,nx - nrepl*n)
  nrep <- rep(nrepl,n)
  nrep[nplus] <- nrepl+1
  x[order(x)] <- rep(seq.int(n),nrep)
  x
}

