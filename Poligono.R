# Calculo de superficies 
install.packages("PBSmapping")
install.packages("polyclip")

library(PBSmapping)
library(polyclip)


#-----------------------------------------------

#1.PREPARACION DE DATOS
#---------------------------------------------------------------------------------------------
#Establece el directorio de trabajo requerido
setwd("C:/Users/mauricio.mardones/Documents/GINGAM/Machas/Eval_dir")
dir()
Areas=c()

#Leer Los datos
Puntos= read.table(file="DensMacha.csv", header=TRUE, sep=";")
View(Puntos)
datos= Puntos[Puntos$Playa=="Godoy",]
head(datos)

#Convertir a formato 
datos.T<-data.frame(PID=rep(1,length(datos[,1])),POS=seq(1,length(datos[,1])),CONTEO=datos$Total,X=datos$X,Y=datos$Y)
head(datos)

#POLIGONO AREA TOTAL
# Obtiene los puntos limites del area
source("generapolyJC.R")#llama a función
nn=((max(datos.T$Y)-min(datos.T$Y))/100);nn
polig   <- genera.polyJC(datos.T$X,datos.T$Y,100,10)
convexH <- data.frame(X=polig[,1],Y=polig[,2])

# Obtains buffered convex hull
poly <- list(x=convexH$X, y=convexH$Y)
aux <- polyoffset(poly, 100, jointype='round')
poly.T <- data.frame(cbind(X=(aux[[1]]$x), Y=aux[[1]]$y))

#Crea el poligono
aux.1 <- rep(unique(datos.T$PID),length(poly.T[,1])); aux.2 <- 1:length(poly.T[,1])
poly.T <- cbind(PID=aux.1, POS=aux.2, X=poly.T[,1],Y=poly.T[,2])
poly.T <- as.PolySet(poly.T); class(poly.T)
poly.T <- fixBound(poly.T, tol=0.000001)
poly.T <- closePolys(poly.T)
plot(poly.T$X ,poly.T$Y, type='l', xlab='Easting (km)', ylab='Northing (km)',main="")
points(datos.T$X,datos.T$Y,cex=.6,pch=20)
write.csv(poly.T, file="P.playa.csv", row.names=FALSE)

#Caluculo de areas
Area.Parche<-calcArea(poly.T)
Areas=rbind(Areas,Area.Parche)

#Guardar los datos
write.csv(Areas,"playaX.CSV", row.names=FALSE)
