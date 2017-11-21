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
#Listado de archivos disponibles
HL01pos<-read.geodata("Corcovado.txt")
View(HL01pos)
#Listado de archivos disponibles
list.files()

View(HL01pos)
dup.coords(HL01pos)
#LC1 <- jitterDupCoords(HL01pos,max=0.001)
names(HL01pos)
plot(HL01pos)
hist(HL01pos$data)


# Los datos no se distribuyen normales y se requiere una #transformacion 
#de la familia Box-Cox lo que es controlado por el parametro "lamda"
#se busca el mejor valor de lambda para la transformacion de los datos
b=boxcox((data) ~ V1 + V2, data = HL01pos, lambda = seq(-5, 5, by=1/300))
lambda=b$x
lik=b$y
bc=cbind(lambda,lik)
bc_ord=bc[order(-lik),]
Lambda_optimo=bc_ord[1,1]
Lambda_optimo

#   lambda 
#0.6466667
plot(HL01pos,lambda=Lambda_optimo)
plot(HL01pos)
summary(HL01pos)
summary(HL01pos,lambda=Lambda_optimo)


#EMPIRICAL VARIOGRAMS
#Variograma de Nubes y variogramas experimentales
cloud1HL01pos <- variog(HL01pos,option="cloud",max.dist=15000)
cloud2HL01pos <- variog(HL01pos,option="cloud",estimator.type="modulus",max.dist=15000)
cloud3HL01pos <- variog(HL01pos,option="cloud",lambda=Lambda_optimo,max.dist=15000)
cloud4HL01pos <- variog(HL01pos,option="cloud",lambda=Lambda_optimo,estimator.type="modulus",max.dist=15000)

par(mfrow = c(2,2))
plot(cloud1HL01pos)
plot(cloud2HL01pos)
plot(cloud3HL01pos)
plot(cloud4HL01pos)


bin1HL01pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=0.65))

bin2HL01pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=0.65),estimator.type = "modulus")

bin1HLlbdNout <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=0.65),lambda=Lambda_optimo)

bin2HLlbdNout <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=0.65),estimator.type = "modulus",lambda=Lambda_optimo)


par(mfrow = c(2,2))
plot(bin1HL01pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin2HL01pos, main = "Modulus estimator",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin1HLlbdNout, main = "Classical estimator lambda",xlab = "Distancia (m)", ylab = "Semivarianza")
plot(bin2HLlbdNout, main = "Modulus estimator lambda",xlab = "Distancia (m)", ylab = "Semivarianza")


# Se busca mejor distancia entre puntos  
bin1HLpos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=15.00),estimator.type = "classical",lambda=Lambda_optimo)

bin1HL2pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=20.00),estimator.type = "classical",lambda=Lambda_optimo)

bin2HLpos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=25.00),estimator.type = "classical",lambda=Lambda_optimo)

bin2HL2pos <- variog(HL01pos, uvec=seq(summary(HL01pos)$distances.summary[["min"]],7500.00,l=30.00),estimator.type = "classical",lambda=Lambda_optimo)

par(mfrow=c(2,2))
plot(bin1HLpos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit<-variofit(bin1HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(HLvario.fit,lwd=2,col="black")
summary(HLvario.fit)

plot(bin1HL2pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit2<-variofit(bin1HL2pos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(HLvario.fit2,lwd=2,col="black")
summary(HLvario.fit2)

plot(bin2HLpos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit3<-variofit(bin2HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(HLvario.fit3,lwd=2,col="black")
summary(HLvario.fit3)                            

plot(bin2HL2pos, main = "Classical estimator",xlab = "Distancia (m)", ylab = "Semivarianza") 
HLvario.fit4<-variofit(bin2HL2pos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(HLvario.fit4,lwd=2,col="black")
summary(HLvario.fit)


# Se ajustan varios modelos al mejor variograma experimental 
    
bin1HLpos  
    
par(mfrow=c(2,2))
plot(bin1HLpos,main="Matern",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit1out<-variofit(bin1HLpos,cov.model="matern",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit1out,lwd=2,col="black")
#legend(x=20, y=3500000, lwd=2, col="black",legend="matern")

plot(bin1HLpos,main="Esférico",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit2out<-variofit(bin1HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit2out,lwd=2,col="black")
#legend(x=70, y=10000, lwd=2, col="red",legend="spherical")

plot(bin1HLpos,main="Gausiano",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit3out<-variofit(bin1HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit3out,lwd=2,col="black")
#legend(x=70, y=10000, lwd=2, col="red",legend="gaussian")

plot(bin1HLpos,main="Exponencial",xlab = "Distancia (m)", ylab = "Semivarianza")
bin1HLvario.fit4out<-variofit(bin1HLpos,cov.model="exponential",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit4out,lwd=2,col="black")
#legend(x=60, y=10000, lwd=2, col="red",legend="exponential")

summary(bin1HLvario.fit1out)
summary(bin1HLvario.fit2out)
summary(bin1HLvario.fit3out)
summary(bin1HLvario.fit4out)

####################################################################################


#Variograma experimental con modelos en la misma figura

plot(bin1HLpos,xlab = "Distancia (m)", ylab = "Semivarianza")

bin1HLvario.fit2out<-variofit(bin1HLpos,cov.model="spherical",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit2out,lwd=3,lty=1,col="grey")

bin1HLvario.fit3out<-variofit(bin1HLpos,cov.model="gaussian",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit3out,lwd=3,lty=2,col="black")

bin1HLvario.fit4out<-variofit(bin1HLpos,cov.model="exponential",fix.nugget=FALSE,max.dist=7500)
lines(bin1HLvario.fit4out,lwd=3,lty=1,col="black")

legend(x=1500, y=30, lty=c(1,2,1), lwd=c(3,3,3), col=c("grey","black", "black"), legend=c("Esférico", "Gausiano", "Exponencial"))





#FIT BY ML                                       
#ini.cov.pars("sill/sigma","rango")#sigma o su equivalente sill      HL01pos  ,lambda=0.025

HL01varioMLmatout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="matern",ini=c(8.3,5.834706),nugget = 727)

# 0.1652176    0.1548648 1956.8899815 
HL01varioMLspheout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="spherical",ini=c(0.1548648,1956.8899815),nugget = 0.1652176)

HL01varioMLgausout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="gaus",ini=c(1.6,354200.4),nugget=1.5)

HL01varioMLexpout<-likfit(HL01pos,lambda=Lambda_optimo,cov.model="exponential",ini=c(2.0,570000),nugget =1.2)


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

#Cross-validation model Matern 
WLS.xval.matHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit2out)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.matHL01pos)

#Cross-validation model spherical 
WLS.xval.spheHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit3out)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.spheHL01pos)

#Cross-validation model gaussian 
WLS.xval.gausHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1LCNvario.fit3out)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.gausHL01pos)

#Cross-validation model exponential 
WLS.xval.expHL01pos<-xvalid(HL01pos,lambda=Lambda_optimo, model=bin1HLvario.fit4out)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(WLS.xval.expHL01pos)

summary(WLS.xval.matHL01pos)
summary(WLS.xval.spheHL01pos)
summary(WLS.xval.gausHL01pos)
summary(WLS.xval.expHL01pos)

##########################################################3
#CROSS-VALIDATION  CON MODELOS ML

summary(HL01varioMLmatout)
summary(HL01varioMLspheout)
summary(HL01varioMLgausout)
summary(HL01varioMLexpout)

#Cross-validation model Matern 
ML.xval.matHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLmatout)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.matHL01)

#Cross-validation model spherical 
ML.xval.spheHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLspheout)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.spheHL01)

#Cross-validation model gaussian 
ML.xval.gausHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLgausout)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.gausHL01)

#Cross-validation model exponential 
ML.xval.expHL01<-xvalid(HL01pos,lambda=Lambda_optimo, model=HL01varioMLexpout)
par(mfcol=c(5, 2), mar=c(3,3,1,0.5), mgp=c(1.5,0.7,0))
plot(ML.xval.expHL01)

summary(ML.xval.matHL01)
summary(ML.xval.spheHL01)
summary(ML.xval.gausHL01)
summary(ML.xval.expHL01)


#######################################################################
#Después de la Validación Cruzada escojo el modelo esferico, ya que###
#tiene el menor valor de CME. Por lo####
#tanto, es el indicado para hacer el kriging###########################
#######################################################################



HL01varioMLspheout$beta.var #varianza de Beta
# [1]   0.004818578

#LCNposnull.lf$info.minimiLCNposnulltion.function



##--STANDARD DEVIATIONS OF ESTIMATES--##

#solve(LCNposnull.lf$info.minimiLCNposnulltion.function$hessian)

##UPDATING PLOTS
x11()
points(HL01pos,pt.divide=c("data.proportional"),lambda=HL01varioMLspheout$lambda,col="gray",yl="Northing (km)",xl="Easting (km)")                    
plot(HL01pos,lambda=HL01varioMLspheout$lambda)


#=====================================================================================
#
# 2.	Poligono para la estimacion

#=====================================================================================

datos<-data.frame(PID=rep(1,length(HL01pos$coords[,1])),POS=seq(1,length(HL01pos$coords[,1]),1),
                  cuadrante=HL01pos$data,E=HL01pos$coords[,1],N=HL01pos$coords[,2])
head(datos)

install.packages("polyclip")
library(polyclip)
source("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir/generapolyJC.R") #llama a función que está en la carpeta Geostat
polig   <- genera.polyJC(datos$E,datos$N,800,100)

Puntos <- locator(type="o", cex=0.5)

Polig_1 <- data.frame(X=c(Puntos$x, Puntos$x[1]), Y=c(Puntos$y, Puntos$y[1]))

x11()
plot(Polig_1 ,type="l")
points(datos$E,datos$N)

write.table(Polig_1, file='H01pospoly1.csv', col.names=TRUE, sep=',')

#############################

#Datos polygono de density positivas 

H01pospoly1<-read.csv("H01pospoly1.csv")
x11()
plot(H01pospoly1,lwd=2,type="s", xlab="Easting (km)", ylab="Northing (km)")

#calcular area
library(PBSmapping)

Areaefecth<-calcArea(polig)
###########################
# MODELOS AJUSTADOS POR ML

HL01varioMLspheout


##--MODELO ESFERICO SPATIAL PREDICTION--## "el kriging se hace solo para comparar y no para estimar"    
## POLIGONO 1
#--en geostatistic verosimilitudista no es neceLCNposnullrio hacer el kriging, solo se compara estimado y observado"--##
                                #Longitud        Latitud 
HL01pos.pred.grid1<-expand.grid(seq(625000,660000,l=500),seq(5875000, 5910000,l=500)) 
HL01pos.krig1<-krige.conv(HL01pos,loc=HL01pos.pred.grid1,krige=krige.control(obj.m=HL01varioMLspheout),borders=H01pospoly1)

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
HL01pos.sphe.BT.mean1<-mean(BCtransform(rnorm(5000,mean=HL01varioMLspheout$beta,sd=sqrt(HL01varioMLspheout$beta.var)),lambda=HL01varioMLspheout$lambda,inverse=TRUE)$data)
HL01pos.sphe.BT.mean1
#[1]  284.5089


HL01pos.sphe.BT.sd1<-sd(BCtransform(rnorm(5000,mean=HL01varioMLspheout$beta,sd=sqrt(HL01varioMLspheout$beta.var)),lambda=HL01varioMLspheout$lambda,inverse=TRUE)$data)
HL01pos.sphe.BT.sd1
#[1]  176.0032

HL01pos.exp.BT.CV1<-HL01pos.sphe.BT.sd1/HL01pos.sphe.BT.mean1
HL01pos.exp.BT.CV1
#[1] 0.6186211

# Coeficiente de variacion = Desv est/ media
HL01pos.sphe.BT.CV1<-(HL01pos.sphe.BT.sd1/beta.substitute1)*100
HL01pos.sphe.BT.CV1
#  22.35209

HL01pos.krig.sd.BT1<-HL01pos.krig.mean1*HL01pos.sphe.BT.CV1
HL01pos.krig.sd.BT1
#[1] 17600.32



