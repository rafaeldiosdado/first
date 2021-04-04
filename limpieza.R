# Instalación de Paquetes
install.packages(c("VIM","DEoptimR","minqa","nloptr","DMwR", "mvoutlier","TTR","caTools",
                   "AppliedPredictiveModeling","caret"),
                 dependencies = c("Depends", "Suggests"))

#########################################################
#  Datos Faltantes                                      #
#########################################################

load("censusn.rda")

#Para ver que columnas tienen valores perdidos
which(colSums(is.na(censusn))!=0)

#Para ver que filas tienen valores perdidos
rmiss=which(rowSums(is.na(censusn))!=0,arr.ind=T)

#Para ver el porcentaje de filas con valores perdidos
length(rmiss)*100/dim(censusn)[1]

#Para ver el porcentaje de valores perdidos en las columnas
colmiss=c(2,6,13)
per.miss.col=100*colSums(is.na(censusn[,colmiss]))/dim(censusn)[1]
per.miss.col

library(VIM)
a=aggr(censusn,numbers=T)
a
summary(a)

#Ejemplo 2
data(tao)
b<-aggr(tao)
b
marginplot(tao[,c("Air.Temp", "Humidity")])

matrixplot(censusn)

# Eliminar casos
census.cl=na.omit(censusn)

# Imputación

library(DMwR)
census<-censusn
for(h in c(2,4,5,6,7,8,9,13,14)){
  census[,h]<-as.factor(census[,h])
}
library(DMwR)
census.c<-centralImputation(census)
census.d<-initialise(census,method="median")

# K-Vecinos más cercanos
census.k<-knnImputation(census)

# IRMI
imputed.tao <- irmi(tao)
summary(imputed.tao)


#########################################################
#  Outliers                                             #
#########################################################

# Univariados
bupa <- read.csv("bupa.txt")
zbupa=cbind(scale(bupa[,-7]),bupa[,7])
zbupa1=zbupa[,1]
rownames(bupa[abs(zbupa1)>2,])
outliers=boxplot(bupa$V1,plot=F)$out
nout=as.character(outliers)
boxplot(bupa$V1,col="blue")
for(i in 1:length(outliers))
{
  text(outliers[i],as.character(which(bupa$V1==outliers[i])),
       cex=.8,pos=4)
}

# Multivariados

# Estimadores robustos
library(mvoutlier)
aq.plot(bupa[bupa$V7==1,1:6],alpha=0.01)

# Densidad Local
library(DMwR)
bupa1=bupa[bupa$V7==1,1:6]
indice=as.numeric(rownames(bupa1))
lof=lofactor(bupa1,10)
lof
indice [order(lof,decreasing=T)][1:10]

# Otros metodos
library(mvoutlier)

# PCOut Method
bupa1=bupa[bupa$V7==1,1:6]
indice=as.numeric(rownames(bupa1))
outlier=pcout(bupa1)
outlier
indice [order(outlier$wfinal,decreasing=F)][1:10]

# Sign Method
outlier1=sign2(bupa1)
outlier1
indice [order(outlier1$x.dist,decreasing=T)][1:10]

# Clusters
library(cluster)
bupa1=bupa[bupa[,7]==1,1:6]
pambupa1=pam(bupa1,20,stand=T)
pambupa1$clusinfo
bupa1[pambupa1$clustering==19,]


#########################################################
#  Transformación                                       #
#########################################################

bupa<-read.table("datos/bupa.txt",header=T,sep=",") 
summary(bupa)

#Z-score
library(reshape)
zbupa<-rescaler(x=bupa[,-7],type="sd")
zbupa<-cbind(zbupa,bupa[,7] )
summary(zbupa)

#Usando la función scale
zbupa<-cbind(scale(bupa[,-7]),bupa[,7])
summary(zbupa)

#Min-Max

# Usando la librería dprep
#source(file="scripts/dprep.R")
#mmbupa<-mmnorm(bupa,minval=0,maxval=1 )[,-7]
#summary(mmbupa)

library(DMwR)
mmbupa=bupa
# for(h in 1:6){
#   mmbupa[,h]=ReScaling(bupa[,h],0,1)
# }
mmbupa[,-7]<-sapply(bupa[,-7],FUN=ReScaling, t.mn=0, t.mx=1)
summary(mmbupa)

#Escalamiento decimal

# Usando dprep
dsbupa<-decscale(bupa)[,-7]
dsbupa<-cbind(dsbupa,bupa[,7])
summary(dsbupa)


#Sigmoidal
sigbupa<-bupa
sigbupa[,-7]<-signorm(bupa[,-7])
summary(sigbupa)

#Haciendo plots para ver el efecto de la normalizacion
par(mfrow=c(1,2))
plot(sort(bupa$V1))
plot(sort(sigbupa$V1))

#Softmax
library(DMwR)
softbupa<-bupa
softbupa[,-7]<-SoftMax(bupa[,-7],lambda=2*pi)
summary(softbupa)


#Gráfico de comparación
par(mfrow=c(2,3))
boxplot(bupa[,1:6],main="bupa")
boxplot(zbupa[,1:6],main="znorm bupa")
boxplot(mmbupa[,1:6],main="min-max bupa")
boxplot(dsbupa[,1:6],main="dec scale bupa")
boxplot(sigbupa[,1:6],main="signorm bupa")
boxplot(softbupa[,1:6],main="softmax bupa")


################################################################################
### Caso de Estudio: Segmentación celular en Screening de alto contenido

library(AppliedPredictiveModeling)
data(segmentationOriginal)

## Retener el conjunto original de entrenamiento 
segTrain <- subset(segmentationOriginal, Case == "Train")

## Remover las tres primeras columnas (identificadores)
segTrainX <- segTrain[, -(1:3)]
segTrainClass <- segTrain$Class

################################################################################
### Transformación de Predictores Individuales

## La columna VarIntenCh3 mide la desviación estándar de la intensidad
## de los pixels en los filamentos de actinina

max(segTrainX$VarIntenCh3)/min(segTrainX$VarIntenCh3)

library(e1071)
skewness(segTrainX$VarIntenCh3)

library(caret)

## Se usa la función preProcess de la librería caret para transformar debido al sesgo
segPP <- preProcess(segTrainX, method = "BoxCox")

## Aplicación de la transformación
segTrainTrans <- predict(segPP, segTrainX)

## Resultados para un predictor 
segPP$bc$VarIntenCh3

par(mfrow=c(1,1))
histogram(~segTrainX$VarIntenCh3,
          xlab = "Unidades Originales",
          type = "count")
histogram(~log(segTrainX$VarIntenCh3),
          xlab = "Unidades en Logaritmos",
          ylab = " ",
          type = "count")

segPP$bc$PerimCh1

histogram(~segTrainX$PerimCh1,
          xlab = "Unidades Originales",
          type = "count")

histogram(~segTrainTrans$PerimCh1,
          xlab = "Data Transformada",
          ylab = " ",
          type = "count")



