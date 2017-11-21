#MODIFICANDO LA LINEA DONDE ESTÁ EL MODELO ESFERICO,
#ELIMINANDO KAPPA, COLOCANDO EL VALOR DEL NUGGET, RANGO
#ESTIMADO CON LAS DENSIDADES POSITIVAS
#LLAMANDO EL POLYGONO 
#PROBANDO CON GRILLA DE CELDAS DE 1 (mejor) ; 0.5 Y 2 KM^2


#--Geostadística verosimilitudista teniendo en cuenta valores positivos y nulos de densidad--
#--Estimation of the effective area occupied by a stock of shrimp (HL01--
#--with survey data and a binomial geostatistical model--

library(geoRglm) #<- for likelihood maximization of the binomial spatial model
library(coda) #<- for Markov Chain diagnostics
library(lattice)
library(geoR)
library(sp)
library(MASS)

setwd("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir")
#ahora con todos los datos, 0 y 1
setwd(wd)
#Listado de archivos disponibles
list.files()
HL01.pos.null1<-read.csv("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir/Corcovado.csv", header = F, sep = ";") #datos solo positivo
View(HL01.pos.null1)
dup.coords(HL01.pos.null1)
HL01.pos.null1
summary(HL01.pos.null1)
x11()
plot(HL01.pos.null1)

#HL01.pos.null1<- subset(HL01.pos.null1, select = c(X, Y, Total))

#Paso 1: Create a matrix with info on each cell, in 7 variables, 
#with all the cells in the grid in a column vector
cell<-100 #spatial cell of 1000 m each side #asignación del tamaño de celdas
min_east<-floor(min(HL01.pos.null1[1]))-cell
max_east<-ceiling(max(HL01.pos.null1[1]))+cell 
min_north<-floor(min(HL01.pos.null1[2]))-cell
max_north<-ceiling(max(HL01.pos.null1[2]))+cell
east<-seq(from = min_east, to =max_east, by= cell) 
north<-seq(from = min_north, to =max_north, by= cell)
grid<-matrix(NA,((length(east)-1)*(length(north)-1)), 7)
    h <- 1
    for (i in 1: (length(east)-1)) {
        for (j in 1: (length(north)-1)) {
            grid[h,]<-c(east[i], east[i+1],north[j], north[j+1],(east[i]+east[i+1])/2,(north[j]+ north[j+1])/2,h)  
            h <- h + 1
        }
    }
#Paso 2:
#Print the grid, transforms to dataframe, and assign names to columns grid
grid<-as.data.frame(grid)
grid
names(grid)<-c("Easting_ini", "Easting_fin", "Northing_ini", "Northing_fin","Easting_avg","Northing_avg","ID")
grid[1:10,]
x11()
plot(grid$Easting_avg,grid$Northing_avg, col="blue", pch=".", xlab="Easting (km)", ylab="Northing (km)")
points(HL01.pos.null1[,1],HL01.pos.null1[,2])

#Paso 3:
#Runs through the grid building a table with the counts (number of observations 
#in the cell, for example number of tows) and the number of positive observations 
#(Density higher than zero) - May take a while ...
    frec<-mat.or.vec(length(grid[,1]),2)
    for (i in 1: length(grid[,1])){
        for (j in 1: length(HL01.pos.null1[,1])){         
            if(HL01.pos.null1[j,1] > grid[i,1] & HL01.pos.null1[j,1] <= grid[i,2]){ 
                if (HL01.pos.null1[j,2] > grid[i,3] & HL01.pos.null1[j,2] <= grid[i,4]){
                    frec[i,1] <-frec[i,1] + 1
                    if (HL01.pos.null1[j,3] > 0) frec[i,2] <-frec[i,2] + 1
                 }           
             }               
        }
    }

#Paso 4:
#Create matrix for geostatistical analysis
HL01.pos.null1bin<-cbind(grid[,5],grid[,6], frec[,1],frec[,2])
HL01.pos.null1bin
plot(HL01.pos.null1bin)
#Paso 5:
#Validation of counts and of positives (to compare later with
#the report from summary)
length(HL01.pos.null1[,1])
#[1] 660 celdas totales en la grilla

sum(HL01.pos.null1bin[,4])
#[1]  477 celdas positivas

#Paso 6:
#Elliminate empty cells (spatial cells that were not posited during the survey)
eval<-HL01.pos.null1bin[,3]>0
HL01.pos.null1bin<-HL01.pos.null1bin[eval,]
HL01.pos.null1bin

#Paso 7:
#Finally, create geodata object
HL01.pos.null1bin<-as.geodata(HL01.pos.null1bin,coords.col= 1:2, units.m.col= 3,data.col=4)
HL01.pos.null1bin
summary(HL01.pos.null1bin)
 
#Number of data points: 103 

#Coordinates summary
#     Coord1  Coord2
#min -490665 5928327
#max 1109335 6528327

#Distance summary
#    min     max 
# 100000 1676305 

#Data summary
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   3.000   5.000   4.631   6.000   9.000 

#Offset variable summary
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   5.000   6.000   6.408   8.000  10.000 



#--END OF DATA PREPOCESSING--
#---------------------------------#
#---------------------------------#

##--Create Borders--## 
HL01poly1<-read.csv("H01pospoly1.csv")
plot(HL01poly1,lwd=2,type="s", xlab="Easting (km)", ylab="Northing (km)")
#--------------------------#-----------------------#

#START OF GEOSTATISTICAL ANALYSIS

# number of attempts o intentos at observing fish stock
HL01.pos.null1bin.n.trials<-sum(HL01.pos.null1bin$units.m) 
HL01.pos.null1bin.n.trials
#[1] 660 #Compare with validation of counts and of positives
HL01.pos.null1bin.n.succ<-sum(HL01.pos.null1bin$data) # number of successes o éxitos
HL01.pos.null1bin.n.succ
#[1] 477 #Compare with validation of counts and of positives
HL01.pos.null1bin.naive<- HL01.pos.null1bin.n.succ/HL01.pos.null1bin.n.trials # probability of success naive non-spatial
HL01.pos.null1bin.naive
#[1]   0.7227273

#Exploratory plot to see the data in spatial-binomial format
#This plot is useful to check that the chosen size of cell 
#is adequate. The cell must be small but not so small that most cells contain only
#one observation in the Trials panel.

#postscript("",width=10,height=10,horizontal=FALSE,onefile=TRUE,paper="special")
par(mfrow=c(1,2),mar=c(2,2,2,1),oma=c(2,2,0,0),bg="white",fg="black")
plot(HL01.pos.null1bin$coords[,1],HL01.pos.null1bin$coords[,2],type="n",xlab="",ylab="",main="")#plot intentos survey_HL01
text(HL01.pos.null1bin$coords[,1],HL01.pos.null1bin$coords[,2],format(HL01.pos.null1bin$units.m),cex=0.6)
rect(HL01.pos.null1bin$coords[,1]-cell/2,HL01.pos.null1bin$coords[,2]-cell/2,HL01.pos.null1bin$coords[,1]+cell/2,HL01.pos.null1bin$coords[,2]+cell/2)
mtext(text="Este (m)",side=1,outer=TRUE)
mtext(text="Norte (m)",side=2,outer=TRUE)
plot(HL01.pos.null1bin$coords[,1],HL01.pos.null1bin$coords[,2],type="n",xlab="",ylab="",main="")#plot éxitos survey_HL01
text(HL01.pos.null1bin$coords[,1],HL01.pos.null1bin$coords[,2],format(HL01.pos.null1bin$data),cex=0.6)
rect(HL01.pos.null1bin$coords[,1]-cell/2,HL01.pos.null1bin$coords[,2]-cell/2,HL01.pos.null1bin$coords[,1]+cell/2,HL01.pos.null1bin$coords[,2]+cell/2)
#mtext(text="Easting (km)",side=1,outer=TRUE)
#mtext(text="Northing (km)",side=2,outer=TRUE)

#al final plot comparativo de intentos y éxitos
#dev.off()

#--Start exploratory analysis for building the Markov Chain (Cadena de Markov)--#
#--This tuning phase shall end with an acceptance rate around--#
#--0.6, si es menor hay que hacer nuevamente; by changing parameters 
#S.scale and thin in mcmcm.control--# 2.1 sill, 6.6 rango de ajuste likelihood                                        
# cov.pars=c(32.773,46.283) nugget=49.058
HL01.pos.null1.spmod<-list(cov.pars=c(6.287,808142.3),beta=1.0,cov.model="spherical",nugget=2.947,family="binomial",link="logit")
HL01.pos.null1.mcmc<-mcmc.control(S.scale=0.001,thin=8) # escala de las iteracciones y thin el numero de iteraciones 8 es 8000
HL01.pos.null1.tune<-glsm.mcmc(HL01.pos.null1bin,model=HL01.pos.null1.spmod,mcmc.input=HL01.pos.null1.mcmc)

# [...]  After a few attempts we are HL01tisfied with S.scale=0.01 and thin=8 which
# leads to an acceptance rate of 95.0%
 

#Now check the mixing with coda
HL01.pos.null1.tune.coda<-create.mcmc.coda(HL01.pos.null1.tune,mcmc.input=HL01.pos.null1.mcmc)
x11()
autocorr.plot(HL01.pos.null1.tune.coda,ask=TRUE)

##--La autocorrelación tiende a cero en todas las celdas--##

#Run likelihood maximiHL01tion
HL01.pos.null1.pre.lf<-prepare.likfit.glsm(HL01.pos.null1.tune)
HL01.pos.null1.lf<-likfit.glsm(HL01.pos.null1.pre.lf,ini.phi=808142.3,hessian=FALSE,nugget.rel=2.947,fix.nugget.rel=FALSE,cov.model="spherical")
#--------------------#-------------------------#
#--------------------#-------------------------#    

#-------------------#----------------------#

summary(HL01.pos.null1.lf)               


HL01.pos.null1.lf

 

#AIC para modelo spatial
(HL01.pos.null1.AICsp<--2*HL01.pos.null1.lf$loglik+2*HL01.pos.null1.lf$npars)
#[1] -103.9906

#NULL model without spatial effect
HL01.pos.null1.lf0<-likfit.glsm(mcmc.obj = HL01.pos.null1.pre.lf, cov.model = "pure.nugget")

summary(HL01.pos.null1.lf0)



#---Akaike Information Criterion(AIC) 
(HL01.pos.null1.AICnonsp<--2*HL01.pos.null1.lf0$loglik+2*HL01.pos.null1.lf0$npars)
#  12.1184

#AIC difference no spatial model and spatial model
HL01.pos.null1.AICsp-HL01.pos.null1.AICnonsp
#[1]-116.109 
 
#---Back transformation of probability of success from the linear link---
#---This estimate is the probability of observing the stock, and it also--- ---
#---is the proportion of the total area covered in the survey that is--- 
#---actually occupied by the stock---


#Antilogit  #P: Probabilidad de observación del stock 
(HL01.pos.null1.lf.p<-exp(HL01.pos.null1.lf$beta)/(1+exp(HL01.pos.null1.lf$beta)))
#[1]  0.5804363

##---Survey area ()---##
#producto del No. celdas y tamaño de las celdas#
(HL01.pos.null1.area<-summary(HL01.pos.null1bin)$n*cell^2)
#[1]   8.9e+11 #unidad del área (m^2)   

##--Effective occupied area of the stock --## 

(HL01.pos.null1.eff.area<-HL01.pos.null1.area*HL01.pos.null1.lf.p)
#[1]  597849338810  #unidad del área del stock de shrimp(m^2)   

#--I haven't found a way to obtain the estimation variance of Sa.pos.null.lf.p--#
#--So I use the mean of the kriging prediction variances as a conservative proxy--#

##-----Kriging-----##
##---Update MCMC (Cadena de Markov)---#
HL01.pos.null1.spmod<-list(cov.pars=c(HL01.pos.null1.lf$cov.pars[1],HL01.pos.null1.lf$cov.pars[2]),beta=HL01.pos.null1.lf$beta,cov.model=HL01.pos.null1.lf$cov.model,nugget=0.5,family="binomial",link="logit")

HL01.pos.null1.mcmc<-mcmc.control(S.scale=0.020,thin=8)
HL01.pos.null1.tune<-glsm.mcmc(HL01.pos.null1bin,model=HL01.pos.null1.spmod,mcmc.input=HL01.pos.null1.mcmc)

#Create a grid to do the kriging at its nodes inside the polygon

HL01.pos.null1pro.pred.grid<-expand.grid(seq(-600000,1200000,l=200),seq(5700000,6800000,l=200))
fix(HL01.pos.null1pro.pred.grid)
#Do the kriging with the update model inside the polygon - Might take a while ...
HL01.pos.null1bin.krig<-glsm.krige(HL01.pos.null1.tune,loc=HL01.pos.null1pro.pred.grid,borders=HL01poly1)
fix(HL01.pos.null1bin.krig)
#glsm.krige: Prediction for a generalised linear spatial model

#Archivo de salida de resultados de interpolación
#Krige.out <- data.frame(rep(1,dim(HL01.pos.null1pro.pred.grid)[1]), HL01.pos.null1pro.pred.grid$x, HL01.pos.null1pro.pred.grid$y, HL01.pos.null1bin.krig$predict, HL01.pos.null1bin.krig$krige.var)
#names(Krige.out) <- c("Foco1","x","y","Dens","Var")
#fix(Krige.out)
#summary(Krige.out)
#hist(Krige.out$Dens)
 
#*******************************************************************************************
#write.table(Krige.out, "1.csv", sep=";", dec=".", row.names=TRUE) #exportar archivo de datos del kriging



image(HL01.pos.null1bin.krig,col=gray(seq(1,0.1,l=40)),xlab="Easting (m)",ylab="Northing (m)")
legend.krige(val=HL01.pos.null1bin.krig$pred,vert=TRUE,x.leg=c(0,100000),y.leg=c(5500000,5800000),col=gray(seq(1,0.1,l=40)))


#colors gradients
rbg.palette<-colorRampPalette(c("white","yellow","orange","red"),space="rgb")
x11()
par(mar=c(0,15,1,0),oma=c(3,0,0,4))
image(HL01.pos.null1bin.krig,col=rbg.palette(50),xlab="Easting (km)",ylab="Northing (km)")
legend.krige(val=HL01.pos.null1bin.krig$pred,vert=TRUE,x.leg=c(750,752),y.leg=c(200,215),col=rbg.palette(50))

#dev.off()

##---Assessment Area---##
mean(HL01.pos.null1bin.krig$pred)
#[1] 0.7440628


(HL01.pos.null1.eff.area.krig<-HL01.pos.null1.area*mean(HL01.pos.null1bin.krig$pred))
#[1] 766384732885

mean(sqrt(HL01.pos.null1bin.krig$krige.var))
#[1]   0.1259757

HL01.pos.null1bin.area<-summary(HL01.pos.null1bin)$n*cell^2
HL01.pos.null1bin.area
#[1] 1.03e+12 unidad del área (km^2)   

##---BT: Biomasa relativa total (como no es una variable contínua es adimensional)---##
#Z: Densidad media estimada con densidades positivas (kriging mean)
#P: Probabilidad de observación del stock de shrimp
#A: Área survey
Z<-787.4128 #unidades de la densidad (h/m^2)
P<-0.5895519    # cuadricula 50000x50000
A<- 8.9e+11 #unidades del área survey (m^2) cuadricula 10000x10000 

(BT<-(Z*P*A)) #Biomasa (kg)
# 

(BT.ton<-BT/1000) #Biomasa (ton)
#   3136.123


#[1]  39801454085 #unidad del área del stock (m^2) 
