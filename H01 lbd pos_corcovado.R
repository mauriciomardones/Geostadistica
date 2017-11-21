library(geoR)
library(sp)
library(MASS)
library(gstat)
library(grid)
library(lattice)
library(survival)
library(splines)
library(Formula)


rm(list=ls(all=TRUE))
#READ DATA
#START
setwd("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir")
#ahora con todos los datos, 0 y 1
dir()
datos<-read.csv("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir/DensMacha.csv", sep = ";")
playa<-c("Blanca","Corcovado","Godoy","MaoToro","Pangal","Pichicuyén","Pucaihuén","Tajamar")

for(i in 1:length(playa)){
dat<-datos[datos$Playa==playa[i]&datos$Total>0,,]
write.table(data.frame(V1=dat$X,V2=dat$Y,data=dat$Total),file=paste("Playa",playa[i],"pos.txt",sep=""),sep=" ",row.names = F,
            col.names = F)}

#Listado de archivos disponibles
HL01pos<-read.geodata(paste("Playa",playa[2],"pos.txt",sep="")) #playa 2 es corcovado
#Listado de archivos disponibles
list.files()

#View(HL01pos)
dup.coords(HL01pos)
#LC1 <- jitterDupCoords(HL01pos,max=0.001)
names(HL01pos)
plot(HL01pos)
hist(HL01pos$data)


# Los datos no se distribuyen normales y se requiere una #transformacion 
#de la familia Box-Cox lo que es controlado por el parametro "lamda"
#se busca el mejor valor de lambda para la transformacion de los datos
b=boxcox((data) ~ V1 + V2, data = HL01pos, lambda = seq(-5, 5, by=1/100))
lambda=b$x
lik=b$y
bc=cbind(lambda,lik)
bc_ord=bc[order(-lik),]
Lambda_optimo=bc_ord[1,1]
Lambda_optimo

#Comentario
#lambda 
#-0.3466667 
x11()
plot(HL01pos,lambda=Lambda_optimo)
summary(HL01pos)
summary(HL01pos,lambda=Lambda_optimo)

#=============================================================================#
#EMPIRICAL VARIOGRAMS                                                         #
#-----------------------------------------------------------------------------#
#variog(geodata, coords = geodata$coords, data = geodata$data,                #
#       uvec = "default", breaks = "default",                                 #
#       trend = "cte", lambda = 1,                                            #
#       option = c("bin", "cloud", "smooth"),                                 #
#       estimator.type = c("classical", "modulus"),                           #
#       nugget.tolerance, max.dist, pairs.min = 2,                            #
#       bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,   #
#       unit.angle = c("radians","degrees"), angles = FALSE, messages, ...)   #
#-----------------------------------------------------------------------------#


####################################################################################

#Variograma de Nubes y variogramas experimentales
# OPTION "cloud"

####################################################################################
summary(HL01pos)$distances.summary[["max"]]

n<-10
distmax<-1500

cloud1HL01pos <- variog(HL01pos,option="cloud",max.dist=distmax)
cloud2HL01pos <- variog(HL01pos,option="cloud",estimator.type="modulus",max.dist=distmax)
cloud3HL01pos <- variog(HL01pos,option="cloud",lambda=Lambda_optimo,max.dist=distmax)
cloud4HL01pos <- variog(HL01pos,option="cloud",lambda=Lambda_optimo,estimator.type="modulus",max.dist=distmax)

X11()
par(mfrow = c(2,2))
plot(cloud1HL01pos)
plot(cloud2HL01pos)
plot(cloud3HL01pos)
plot(cloud4HL01pos)


####################################################################################
# variogramas experimentales
# OPTION "bin"

####################################################################################

bin1HL01pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=n))
bin2HL01pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=n),estimator.type = "modulus")
bin1HLlbdNout <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=n),lambda=Lambda_optimo)
bin2HLlbdNout <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=n),estimator.type = "modulus",lambda=Lambda_optimo)

x11()
par(mfrow = c(2,2))
plot(bin1HL01pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin2HL01pos, main = "Modulus estimator",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin1HLlbdNout, main = "Classical estimator lambda",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin2HLlbdNout, main = "Modulus estimator lambda",xlab = "Distancia (m)", ylab = "Semivarianza")


####################################################################################

# Se busca mejor distancia entre puntos  VARIOGRAMA ESFERICO

####################################################################################

bin1HLpos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=10.00),estimator.type = "classical",lambda=Lambda_optimo)
bin1HL2pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=15.00),estimator.type = "classical",lambda=Lambda_optimo)
bin2HLpos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=20.00),estimator.type = "classical",lambda=Lambda_optimo)
bin2HL2pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=25.00),estimator.type = "classical",lambda=Lambda_optimo)

x11()
par(mfrow=c(2,2))
plot(bin1HLpos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit<-variofit(bin1HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(HLvario.fit,lwd=2,col="black")
summary(HLvario.fit)

method      <-summary(HLvario.fit)$pmethod
covmodel    <-summary(HLvario.fit)$cov.model
sigma       <-summary(HLvario.fit)$spatial.component[1]
phi         <-summary(HLvario.fit)$spatial.component[2]
kappa       <-summary(HLvario.fit)$spatial.component.extra
nugget      <-summary(HLvario.fit)$nugget.component
rango       <-summary(HLvario.fit)$practicalRange
SumofSquares<-summary(HLvario.fit)$sum.of.squares

m1<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)

plot(bin1HL2pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit2<-variofit(bin1HL2pos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(HLvario.fit2,lwd=2,col="black")
summary(HLvario.fit2)

method      <-summary(HLvario.fit2)$pmethod
covmodel    <-summary(HLvario.fit2)$cov.model
sigma       <-summary(HLvario.fit2)$spatial.component[1]
phi         <-summary(HLvario.fit2)$spatial.component[2]
kappa       <-summary(HLvario.fit2)$spatial.component.extra
nugget      <-summary(HLvario.fit2)$nugget.component
rango       <-summary(HLvario.fit2)$practicalRange
SumofSquares<-summary(HLvario.fit2)$sum.of.squares

m2<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)

plot(bin2HLpos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit3<-variofit(bin2HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(HLvario.fit3,lwd=2,col="black")
summary(HLvario.fit3)                            

method      <-summary(HLvario.fit3)$pmethod
covmodel    <-summary(HLvario.fit3)$cov.model
sigma       <-summary(HLvario.fit3)$spatial.component[1]
phi         <-summary(HLvario.fit3)$spatial.component[2]
kappa       <-summary(HLvario.fit3)$spatial.component.extra
nugget      <-summary(HLvario.fit3)$nugget.component
rango       <-summary(HLvario.fit3)$practicalRange
SumofSquares<-summary(HLvario.fit3)$sum.of.squares

m3<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)

plot(bin2HL2pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit4<-variofit(bin2HL2pos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(HLvario.fit4,lwd=2,col="black")
summary(HLvario.fit4)
method      <-summary(HLvario.fit4)$pmethod
covmodel    <-summary(HLvario.fit4)$cov.model
sigma       <-summary(HLvario.fit4)$spatial.component[1]
phi         <-summary(HLvario.fit4)$spatial.component[2]
kappa       <-summary(HLvario.fit4)$spatial.component.extra
nugget      <-summary(HLvario.fit4)$nugget.component
rango       <-summary(HLvario.fit4)$practicalRange
SumofSquares<-summary(HLvario.fit4)$sum.of.squares

m4<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)


data.frame(m1,m2,m3,m4)

# variograma direccional
bin1HLposd <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],distmax,l=10.00),estimator.type = "modulus",direction=pi/4)
x11()
par(mfrow=c(2,2))
plot(bin1HLposd, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvariod.fit<-variofit(bin1HLposd,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(HLvariod.fit,lwd=2,col="black")
summary(HLvariod.fit)

method      <-summary(HLvariod.fit)$pmethod
covmodel    <-summary(HLvariod.fit)$cov.model
sigma       <-summary(HLvariod.fit)$spatial.component[1]
phi         <-summary(HLvariod.fit)$spatial.component[2]
kappa       <-summary(HLvariod.fit)$spatial.component.extra
nugget      <-summary(HLvariod.fit)$nugget.component
rango       <-summary(HLvariod.fit)$practicalRange
SumofSquares<-summary(HLvariod.fit)$sum.of.squares

m1d<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)



####################################################################################

# Se ajustan varios modelos al mejor variograma experimental 

####################################################################################

bin1HLpos  
    
par(mfrow=c(2,2))
plot(bin1HLpos,main="Matern",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit1out<-variofit(bin1HLpos,cov.model="matern",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit1out,lwd=2,col="black")
#legend(x=20, y=3500000, lwd=2, col="black",legend="matern")

summary(bin1HLvario.fit1out)

method      <-summary(bin1HLvario.fit1out)$pmethod
covmodel    <-summary(bin1HLvario.fit1out)$cov.model
sigma       <-summary(bin1HLvario.fit1out)$spatial.component[1]
phi         <-summary(bin1HLvario.fit1out)$spatial.component[2]
kappa       <-summary(bin1HLvario.fit1out)$spatial.component.extra
nugget      <-summary(bin1HLvario.fit1out)$nugget.component
rango       <-summary(bin1HLvario.fit1out)$practicalRange
SumofSquares<-summary(bin1HLvario.fit1out)$sum.of.squares

mb1<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)

plot(bin1HLpos,main="Esférico",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit2out<-variofit(bin1HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit2out,lwd=2,col="black")
#legend(x=70, y=10000, lwd=2, col="red",legend="spherical")

summary(bin1HLvario.fit2out)
method      <-summary(bin1HLvario.fit2out)$pmethod
covmodel    <-summary(bin1HLvario.fit2out)$cov.model
sigma       <-summary(bin1HLvario.fit2out)$spatial.component[1]
phi         <-summary(bin1HLvario.fit2out)$spatial.component[2]
kappa       <-summary(bin1HLvario.fit2out)$spatial.component.extra
nugget      <-summary(bin1HLvario.fit2out)$nugget.component
rango       <-summary(bin1HLvario.fit2out)$practicalRange
SumofSquares<-summary(bin1HLvario.fit2out)$sum.of.squares

mb2<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)

plot(bin1HLpos,main="Gausiano",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit3out<-variofit(bin1HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit3out,lwd=2,col="black")
#legend(x=70, y=10000, lwd=2, col="red",legend="gaussian")

summary(bin1HLvario.fit3out)
method      <-summary(bin1HLvario.fit3out)$pmethod
covmodel    <-summary(bin1HLvario.fit3out)$cov.model
sigma       <-summary(bin1HLvario.fit3out)$spatial.component[1]
phi         <-summary(bin1HLvario.fit3out)$spatial.component[2]
kappa       <-summary(bin1HLvario.fit3out)$spatial.component.extra
nugget      <-summary(bin1HLvario.fit3out)$nugget.component
rango       <-summary(bin1HLvario.fit3out)$practicalRange
SumofSquares<-summary(bin1HLvario.fit3out)$sum.of.squares

mb3<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)


plot(bin1HLpos,main="Exponencial",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit4out<-variofit(bin1HLpos,cov.model="exponential",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit4out,lwd=2,col="black")
#legend(x=60, y=10000, lwd=2, col="red",legend="exponential")
summary(bin1HLvario.fit4out)
method      <-summary(bin1HLvario.fit4out)$pmethod
covmodel    <-summary(bin1HLvario.fit4out)$cov.model
sigma       <-summary(bin1HLvario.fit4out)$spatial.component[1]
phi         <-summary(bin1HLvario.fit4out)$spatial.component[2]
kappa       <-summary(bin1HLvario.fit4out)$spatial.component.extra
nugget      <-summary(bin1HLvario.fit4out)$nugget.component
rango       <-summary(bin1HLvario.fit4out)$practicalRange
SumofSquares<-summary(bin1HLvario.fit4out)$sum.of.squares

mb4<-rbind(method,covmodel,sigma,phi,kappa,nugget,rango,SumofSquares)


data.frame(mb1,mb2,mb3,mb4)

####################################################################################

# Variograma experimental con modelos en la misma figura

####################################################################################
x11()
plot(bin1HLpos,xlab = "Distancia (m)", ylab = "Semivarianza")

bin1HLvario.fit2out<-variofit(bin1HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit2out,lwd=3,lty=1,col="grey")

bin1HLvario.fit3out<-variofit(bin1HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit3out,lwd=3,lty=2,col="black")

bin1HLvario.fit4out<-variofit(bin1HLpos,cov.model="exponential",fix.nugget=FALSE,max.dist=distmax)
lines(bin1HLvario.fit4out,lwd=3,lty=1,col="black")

legend(x=10, y=32, lty=c(1,2,1), lwd=c(3,3,3), col=c("grey","black", "black"), legend=c("Esférico", "Gausiano", "Exponencial"))

####################################################################################

# FIT BY ML    

####################################################################################

#ini.cov.pars("sill/sigma","rango")#sigma o su equivalente sill      HL01pos  ,lambda=0.025

HL01varioMLmatout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="matern",ini=c(9.45143,1085.063),nugget = 6.475415)

HL01varioMLspheout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="spherical",ini=c(8.574355,893.7985),nugget = 7.042699)

HL01varioMLgausout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="gaus",ini=c(7.580788,787.5544),nugget=8.176406)

HL01varioMLexpout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="exponential",ini=c(9.45143,1085.063),nugget =6.475415)


plot(bin1HLpos,xlab = "Distancia (m)", ylab = "Semivarianza")
plot(cloud3HL01pos)
lines(HL01varioMLmatout, col="gray",lwd=2)
lines(HL01varioMLspheout,col="black",lwd=2)
lines(HL01varioMLgausout, col="red",lwd=2)
lines(HL01varioMLexpout, col="blue",lwd=2)

summary(HL01varioMLmatout)
summary(HL01varioMLspheout)
summary(HL01varioMLgausout)
summary(HL01varioMLexpout)


##############################################################

#CROSS-VALIDATION  CON MODELOS WLS

##############################################################
#Cross-validation model Matern 
WLS.xval.matHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit1out)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.matHL01pos)

#Cross-validation model spherical 
WLS.xval.spheHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit2out)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.spheHL01pos)

#Cross-validation model gaussian 
WLS.xval.gausHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit3out)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.gausHL01pos)

#Cross-validation model exponential 
WLS.xval.expHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit4out)

x11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.expHL01pos)

summary(WLS.xval.matHL01pos)
summary(WLS.xval.spheHL01pos)
summary(WLS.xval.gausHL01pos)
summary(WLS.xval.expHL01pos)

##########################################################

#CROSS-VALIDATION  CON MODELOS ML

##########################################################

summary(HL01varioMLmatout)
summary(HL01varioMLspheout)
summary(HL01varioMLgausout)
summary(HL01varioMLexpout)

#Cross-validation model Matern 
ML.xval.matHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLmatout)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.matHL01)

#Cross-validation model spherical 
ML.xval.spheHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLspheout)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.spheHL01)

#Cross-validation model gaussian 
ML.xval.gausHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLgausout)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.gausHL01)

#Cross-validation model exponential 
ML.xval.expHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLexpout)
X11();par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.expHL01)

summary(ML.xval.matHL01)
summary(ML.xval.spheHL01)
summary(ML.xval.gausHL01)
summary(ML.xval.expHL01)


#######################################################################
#Después de la Validación Cruzada escojo el modelo gausiano, ya que###
#tiene el menor valor de CME. Por lo####
#tanto, es el indicado para hacer el kriging###########################
#######################################################################



HL01varioMLgausout$beta.var #varianza de Beta
# [1]   0.004818578

#LCNposnull.lf$info.minimiLCNposnulltion.function

##--STANDARD DEVIATIONS OF ESTIMATES--##

#solve(LCNposnull.lf$info.minimiLCNposnulltion.function$hessian)

##UPDATING PLOTS
x11()
points(HL01pos,pt.divide=c("data.proportional"),lambda=HL01varioMLgausout$lambda,col="gray",yl="Northing (km)",xl="Easting (km)")                    
plot(HL01pos,lambda=HL01varioMLgausout$lambda)


#=====================================================================================
#
# 2.	Poligono para la estimacion

#=====================================================================================

#Leer Los datos


#aca creo solo se leen los positivos

Puntos= read.table(file="DensMacha.csv", header=TRUE, sep=";")
View(Puntos)
dat2<-datos[datos$Playa==playa[2],,]

head(datos)

dat<-data.frame(PID=rep(1,length(dat2$X),1),POS=seq(1,length(dat2$Y),1),
                  cuadrante=dat2$Total,E=dat2$X,N=dat2$Y)
head(dat)

#install.packages("polyclip")
library(polyclip)
library(PBSmapping)
#POLIGONO AREA TOTAL
# Obtiene los puntos limites del area
source("generapolyJC.R")#llama a función
nn=((max(dat$N)-min(dat$N))/100);nn
polig   <- genera.polyJC(dat$E,dat$N,100,10)
convexH <- data.frame(E=polig[,1],N=polig[,2])

# Obtains buffered convex hull
poly <- list(x=convexH$E, y=convexH$N)
aux <- polyoffset(poly, 100, jointype='round')
poly.T <- data.frame(cbind(E=(aux[[1]]$x), N=aux[[1]]$y))

#Crea el poligono
aux.1 <- rep(unique(dat$PID),length(poly.T[,1])); aux.2 <- 1:length(poly.T[,1])
poly.T <- cbind(PID=aux.1, POS=aux.2, X=poly.T[,1],Y=poly.T[,2])
poly.T <- as.PolySet(poly.T); class(poly.T)
poly.T <- fixBound(poly.T, tol=0.00001)
poly.T <- closePolys(poly.T)
X11()
plot(poly.T$X ,poly.T$Y, type='l', xlab='Easting (km)', ylab='Northing (km)',main="")
points(dat$E,dat$N,cex=1,pch=20)

write.csv(poly.T, file="PoliPlayaCorcovado.csv", row.names=FALSE)

#write.table(data.frame(V1=poly.T$X,V2=poly.T$Y),file="PoliPlayaCorcovado.txt",sep=" ",row.names = F,
#            col.names = F)

#Caluculo de areas
Area.Parche<-calcArea(poly.T)

########################################################################################################

# MODELOS AJUSTADOS POR ML

########################################################################################################

HL01varioMLgausout

H01pospoly1<-read.csv("PoliPlayaCorcovado.csv", sep=";")

max(H01pospoly1$x)
min(H01pospoly1$x)
max(H01pospoly1$y)
min(H01pospoly1$y)

########################################################################################################

##--MODELO ESFERICO SPATIAL PREDICTION--## "el kriging se hace solo para comparar y no para estimar"    

########################################################################################################
## POLIGONO 1
#--en geostatistic verosimilitudista no es neceLCNposnullrio hacer el kriging, solo se compara estimado y observado"--##
                                #Longitud        Latitud 
HL01pos.pred.grid1<-expand.grid(seq(min(H01pospoly1$x),max(H01pospoly1$x),l=200),seq(min(H01pospoly1$y), max(H01pospoly1$y),l=200)) 
HL01pos.krig1<-krige.conv(HL01pos,loc=HL01pos.pred.grid1,krige=krige.control(obj.m=HL01varioMLgausout),borders=H01pospoly1)


### HASTA AQUI!!!! error en el poligono, no se si no lo lee bien o que, 
#parece que tiene que ser un poligono de los positivos o total,, revisar!!!!



#postscript(width=5,height=5,horizontal=FALSE,onefile=TRUE,paper="special")
X11()
image(HL01pos.krig1,loc=HL01pos.pred.grid1,col=gray(seq(1,0.1,l=30)),xlab="Easting (m)",ylab="Northing (m)")
legend.krige(val=HL01pos.krig1$pred,vert=TRUE,x.leg=c(0,100000),y.leg=c(630000,5905000),col=gray(seq(1,0.1,l=8)),zlim=c(0,max(HL01pos.krig1$pred)))


##MEAN PROCESS
beta.substitute1<-HL01pos.krig.mean1<-mean(HL01pos.krig1$predict)
#Kriging mean (density mean)# 
beta.substitute1
#[1]  2.747758 #Densidad media

##Backtransform of beta using Montecarlo, to obtain the coefficient of variation
HL01pos.gaus.BT.mean1<-mean(BCtransform(rnorm(5000,mean=HL01varioMLgausout$beta,sd=sqrt(HL01varioMLgausout$beta.var)),lambda=HL01varioMLgausout$lambda,inverse=TRUE)$data)
HL01pos.gaus.BT.mean1
#[1]  284.5089


HL01pos.gaus.BT.sd1<-sd(BCtransform(rnorm(5000,mean=HL01varioMLgausout$beta,sd=sqrt(HL01varioMLgausout$beta.var)),lambda=HL01varioMLgausout$lambda,inverse=TRUE)$data)
HL01pos.gaus.BT.sd1
#[1]  176.0032

HL01pos.gaus.BT.CV1<-HL01pos.gaus.BT.sd1/HL01pos.gaus.BT.mean1
HL01pos.gaus.BT.CV1
#[1] 0.6186211

# Coeficiente de variacion = Desv est/ media
HL01pos.gaus.BT.CV1<-(HL01pos.gaus.BT.sd1/beta.substitute1)*100
HL01pos.gaus.BT.CV1
#  22.35209

HL01pos.krig.sd.BT1<-HL01pos.krig.mean1*HL01pos.gaus.BT.CV1
HL01pos.krig.sd.BT1
#[1] 17600.32



rbind(beta.substitute1,HL01pos.gaus.BT.mean1,HL01pos.gaus.BT.sd1,HL01pos.gaus.BT.CV1,HL01pos.krig.sd.BT1)


