library(RGeostats)
library(geoR)
library(gstat)
library(geospt)
library(geosptdb)
library(sp)
library(sf)
#library(intamap) # Revisar, aparentemente fue descontinuada
library(FNN)
library(xtable)
library(readxl)
library(ggplot2)
library(scatterplot3d)
library(ggplot2)
library(dplyr)
library(tidyr)
library(lattice)
library(gridExtra)
library(beepr)


rm(list = ls())
##############################################
#                    Data                    #
##############################################
#--------------------------------------------#
#               Data Portatil                #
#--------------------------------------------#

setwd("D:/programacion/ejerciciosRGeoestadistica/proyectoFinal")
path <- "../proyectoFinal/data/PluviDBFinal.csv"
pathCundi <- "../proyectoFinal/data/munCundiShp/Municipios_Cundinamarca.shp"
data <- read.csv(path)
save.image(file="../proyectoFinal/.Rdata")
load("../proyectoFinal/.Rdata")


#---------------------------------------------#
#              Data Computador                #
#---------------------------------------------#

setwd("D:/Universidad/7 semestre/Geoestadistica/Ejercicios/ejerciciosRGeoestadistica/proyectoFinal")
path <- "../proyectoFinal/data/PluvidbFinal.csv"
pathCundi <- "../proyectoFinal/data/munCundiShp/Municipios_Cundinamarca.shp"
data <- read.csv(path)
save.image("../proyectoFinal/.Rdata")
load("../proyectoFinal/.Rdata")

##############################################
#             Funciones Externas             #
##############################################
# Calcula los poligonos de Voronoi
source("../proyectoFinal/scriptsR/polVor.R")
# Calcula la media por declustering
source("../proyectoFinal/scriptsR/meanDeclu.R")
# Calcula el variograma Experimental a partir de est.variograms
source("../proyectoFinal/scriptsR/varExp.R")
# Ajusta una lista de modelos por Minimos cuadrados
source("../proyectoFinal/scriptsR/adjustMc.R")
# Ajusta mediante ML y REML
source("../proyectoFinal/scriptsR/adjustML.R")
# Plotea Ajustes Obtenidos por Mínimos Cuadrados
source("../proyectoFinal/scriptsR/plotAdVarMc.R")
# Plotea Ajustes Obtenidos por Máxima Verosimilitud
source("../proyectoFinal/scriptsR/plotAdML.R")
##############################################
#       Manipulación Base de Datos           #
##############################################

# Shape Cundinamarca y sus Municipios
cundi <- st_read(pathCundi)
cundi <- st_transform(cundi, crs = 3116)
cundi.sp <- as(cundi, "Spatial")
proj4string(cundi.sp) <- CRS("+init=epsg:3116")
x11()

# Cambio de Formato datos observados
df <- dplyr::select(data, Nombre, Latitud, Longitud, Municipio, starts_with("X"))
for (col in names(df)[5:ncol(df)]) {
  df[[col]] <- as.numeric(gsub(",", ".", df[[col]]))
}

# Creación del DataFrame a Trabajar
preci <- dplyr::select(df, Nombre, Latitud,Longitud,Municipio,"X2020.01")
preci[preci == 0] <- NA
preci <- na.omit(preci)
preci<- st_as_sf(preci, coords = c("Longitud", "Latitud"), crs = 4326)
preci <- st_transform(preci, crs = 3116)
preci <- preci %>% mutate(Municipio = toupper(Municipio))
preci <- preci %>% rename("z" = "X2020.01")

preci.df <- data.frame()
preci.df <- cbind(st_drop_geometry(preci), st_coordinates(preci))
preci.df <- preci.df[ , !(names(preci.df) %in% "Nombre") ]
preci.df$id <- preci.df$id <- seq_len(nrow(preci.df))
preci.df <- preci.df %>% rename (x=X,y=Y,mpio=Municipio)
preci.df <- preci.df[,sort(names(preci.df))]

# Creación del DataFrame Espacial
preci.sp <- preci.df
coordinates(preci.sp) <- ~x+y
geoPreci <- as.geodata(preci.df, coords.col = 3:4, data.col = 5)
plot(geoPreci)

# Resumen de estadísticas
summary(geoPreci)

# Histograma 
hist(geoPreci$data, freq = FALSE, breaks = 10, xlab = "Precipitación (mm)", ylab = "Frecuencia",
     main="Histograma de la Precipitación")
curve(dnorm(x, mean = mean(geoPreci$data), sd = sd(geoPreci$data)), add = TRUE)

# Q-Q Plot
qqnorm(geoPreci$data, ylab = "Precipitación (mm)", xlab = "Cuantiles teóricos",
       main="QQ-Plot Precipitación")
qqline(geoPreci$data)

# Box-Plot
boxplot(geoPreci$data, main = "Box-Plot Precipitación", notch = FALSE, horizontal = TRUE, xlab = "Precipitación (mm)")

# Pruebas de normalidad
shapiro.test(geoPreci$data) 

# Obteniendo Distancia Máxima
mtxDist <- as.matrix(dist(geoPreci$coords))
rango <-  max(mtxDist)/2

# Verificando Tendencia
modelo1 <- lm(z~x+y,data=preci.df)
summary(modelo1)

# --- Convertir el shapefile de municipios a Spatial y asignar proyección ---
cundiSp <- as(cundi, "Spatial")
proj4string(cundiSp) <- CRS("+init=epsg:3116")

# --- Generar la grilla regular de puntos dentro del shape ---
puntos <- spsample(cundiSp, n = 30000, type = "regular")
gridded(puntos) <- TRUE
proj4string(puntos) <- CRS("+init=epsg:3116")

##############################################
#                    IDW                     #
##############################################

# --- Buscar el mejor valor de q  ---
pVal <- seq(0.25, 10, 0.25)
rmspeVal <- numeric(length(pVal))
nVal <- list(20,10)
pb <- txtProgressBar(min = 0, max = length(pVal), style = 3)
for (i in seq_along(pVal)) {
  p <- pVal[i]
  rmspe <- tryCatch({
    idw.cv(z ~ 1, ~x + y, data = preci.df, nmax = nVal[[1]], nmin = nVal[[2]], p = p)
  }, error = function(e) NA)
  rmspeVal[i] <- rmspe
  setTxtProgressBar(pb, i)
}
close(pb) 

# --- Crear SpatialPointsDataFrame para df1 con coordenadas ---
proj4string(preci.sp) <- CRS("+init=epsg:3116")

# --- Interpolación IDW con el mejor q ---
interIdw <- idw(z ~ 1, preci.sp, newdata = puntos, idp = pVal[which.min(rmspeVal)], nmax = nVal[[1]], nmin = nVal[[2]])
gridded(interIdw) <-  TRUE
# --- Visualización de la interpolación ---
spplot(interIdw["var1.pred"], cuts = 60, scales = list(draw = TRUE), 
       xlab = "Este (m)", ylab = "Norte (m)", main = "Interpolación IDW de Precipitación", 
       auto.key = FALSE)

##############################################
#                    RBF                     #
##############################################

#Multicuadrática 
preciM <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=64,func="M",eta.dmax=2,rho.dmax=2,x0=c(0.1,.01),iter=150)
predRbfM <- rbf(z~1, preci.sp, eta=preciM$Opt$par[1], rho=preciM$Opt$par[2], newdata=puntos, n.neigh=64, func="M")
coordinates(predRbfM) <- c("x", "y")
gridded(predRbfM) <- TRUE
spplot(predRbfM["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Multicuadrático (M)", key.space=list(space="right", cex=0.8))
beep()

# Inversa Multicuadrática
preciIm <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=5,func="IM",eta.dmax=30,rho.dmax=30,x0=c(0.1,.01),iter=500)
predRbfIM <- rbf(z~1, preci.sp, eta=preciIm$Opt$par[1], rho=preciIm$Opt$par[2], newdata=puntos, n.neigh=5, func="IM")
coordinates(predRbfIM)<- c("x", "y")
gridded(predRbfIM)<- TRUE
spplot(predRbfIM["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Inversa Multicuadrática (IM)", key.space=list(space="right", cex=0.8))
beep()

# Spline con Tensión
preciSt <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=32,func="ST",eta.dmax=2,rho.dmax=2,x0=c(0.1,.01),iter=250)
predRbfST <- rbf(z~1, preci.sp, eta=preciSt$Opt$par[1], rho=preciSt$Opt$par[2], newdata=puntos, n.neigh=32, func="ST")
coordinates(predRbfST) <- c("x", "y")
gridded(predRbfST) <- TRUE
spplot(predRbfST["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Spline Con Tensión (ST)", key.space=list(space="right", cex=0.8))
beep()

# Spline de Capa Delgada 
preciTps <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=99,func="TPS", eta.dmax=2,rho.dmax=2,x0=c(0.1,.01),iter=100)
predRbfTPS <- rbf(z~1, preci.sp, eta=preciTps$Opt$par[1], rho=preciTps$Opt$par[2], newdata=puntos, n.neigh=99, func="TPS")
coordinates(predRbfTPS) <- c("x", "y")
gridded(predRbfTPS) <- TRUE
spplot(predRbfTPS["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Spline de Capa Delgada (TPS)", key.space=list(space="right", cex=0.8))
beep()

# Spline Completamente Regularizado
preciCrs <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=48,func="CRS", eta.dmax=2,rho.dmax=2,x0=c(0.1,.01),iter=500)
predRbfCRS <- rbf(z~1, preci.sp, eta=preciCrs$Opt$par[1], rho=preciCrs$Opt$par[2], newdata=puntos, n.neigh=48, func="CRS")
coordinates(predRbfCRS) <- c("x", "y")
gridded(predRbfCRS) <- TRUE
spplot(predRbfCRS["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Spline Completamente Regularizado (CRS)", key.space=list(space="right", cex=0.8))
beep()

# Gaussiano
preciGau <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=5,func="GAU", eta.dmax=10,rho.dmax=10,x0=c(0.1,.01),iter=500)
predRbfGAU <- rbf(z~1, preci.sp, eta=preciGau$Opt$par[1], rho=preciGau$Opt$par[2] , newdata=puntos, n.neigh=5, func="GAU")
coordinates(predRbfGAU) <- c("x", "y")
gridded(predRbfGAU) <- TRUE
spplot(predRbfGAU["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Gaussiano", key.space=list(space="right", cex=0.8))
beep()

# Exponencial
preciExp <- graph.rbf(z~1,preci.sp,eta.opt=TRUE,rho.opt=TRUE,n.neigh=80,func="EXPON", eta.dmax=2,rho.dmax=2,x0=c(0.1,0.1),iter=500)
predRbfEXPON <- rbf(z~1, preci.sp, eta=preciExp$Opt$par[1], rho=preciExp$Opt$par[2], newdata=puntos, n.neigh=80, func="EXPON")
coordinates(predRbfEXPON)   <- c("x", "y")
gridded(predRbfEXPON)   <- TRUE
spplot(predRbfEXPON["var1.pred"], cuts=40, col.regions=bpy.colors(100), main="Exponencial", key.space=list(space="right", cex=0.8))
beep()

##############################################
#                  Kriging                   #
##############################################
#--------------------------------------------#
#              Kriging Simple                #
#--------------------------------------------#

# Obteniendo Normalidad
mdb <- preci.df
mdb.rgdb <- db.create(mdb,ndim=2,autoname=F)
mdb.herm <- anam.fit(mdb.rgdb, name="z", type="gaus") # Name se refiere a la variable a transformar
mdb.hermtrans <- anam.z2y(mdb.rgdb, names="z", anam=mdb.herm)
preci.trans <- as.data.frame(mdb.hermtrans@items)
preci.trans$z <- preci.trans$Gaussian.z
preci.trans <- preci.trans%>% select(-Gaussian.z,-rank)

# Verificación normalidad
plot(density(preci.df$z), main="Distribución df1 Observados")
plot(density(preci.trans$z), main="Distribución df1 Transformados")
shapiro.test(preci.trans$z)
geoPreciTrans <- as.geodata(preci.trans, coords.col=3:4, data.col = 5) 

# Proceso de Declustering
library(spatial)
library(sgeostat)
namesDfs <- c("preci.df", "preci.sp", "geoPreci", "geoPreciTrans", "preci.trans")
res <- meanDeclu(geoPreci$coords,list(preci.df, preci.sp, geoPreci, geoPreciTrans, preci.trans), 
                 knn = 15, namesDfs)
list2env(setNames(res$dfs, namesDfs), envir = .GlobalEnv)
medValues <- res$tabla
print(medValues)

# Verificando anisotropia de los Errores
# estimateAnisotropy(preci.sp, "Err")  -- > Buscar remplazo
var4 <- variog4(geoPreciTrans, data = geoPreciTrans$Err, max.dist=rango)
plot(var4, col = c("red", "blue", "green", "purple"), lwd = 2)

# Calculo  del Error
lags <- 15
preci.pts <- point(preci.trans, x="x", y ="y")
preci.pair <- pair(preci.pts, num.lags=lags,maxdist=rango)
preci.vars <- est.variograms(preci.pts,preci.pair,"Err",trim=0.1)
ylim <- range(preci.vars[, c("robust", "med", "classic", "trimmed.mean")], na.rm = TRUE)
plot(preci.vars$bins,preci.vars$robust,lty=1, col=1,main = "Modelos Semivarianzas Empiricos",xlab="Distancia", ylab="Semivarianza", type="l",ylim=ylim)
lines(preci.vars$bins,preci.vars$med, col=2)
lines(preci.vars$bins,preci.vars$classic, col=3)
lines(preci.vars$bins,preci.vars$trimmed.mean, col=4)
legend("bottomright", c("Robusto", "Mediana", "Clásico","Media Recortada"), col=c(1,2,3,4), lty=c(1,1,1,1))

detach("package:spatial")
detach("package:sgeostat")

# Variograma Experimental con variogram
clasicVar <- variogram(z~1,~x+y,data = preci.trans, beta = medValues$Ponderador[5], cutoff=rango)

# Variograma Experimental para Kriging Simple con est.variograms
ksExp <- varExp(preci.vars,lags,"classic")
hSeq <- seq(0, max(ksExp$dist, na.rm = TRUE), length.out = 15)
varExpKs <- variog(geodata=geoPreciTrans,coords=geoPreciTrans$coords,data=geoPreciTrans$Err, 
                   trend="cte", uvec=hSeq,option="bin",max.dist=max(ksExp$dist), estimator.type = "modulus")

# Ploteo Ambos Variogramas
plot(ksExp$dist, ksExp$gamma, main = expression(paste("Comparación de ", gamma(h), " experimental")),
     xlab = expression("Distancia " ~ h), ylab = expression(gamma(h)), type="l",
     lwd=2, col="black",xlim=c(0,max(c(ksExp$dist, varExpKs$u), na.rm = TRUE)),ylim=c(0,max(c(ksExp$gamma, varExpKs$v), na.rm = TRUE)))
lines(varExpKs$u,varExpKs$v,col="red",type="l",lwd=2)
lines(clasicVar$dist, clasicVar$gamma,col="purple",type="l", lwd="2")
legend("bottomright", legend = c("Variograma con est.variograms", "Variogramas con variog","Clásico con Variogram"),
       col = c("black", "red","purple"), lty = 1, lwd = 2)

# Ajuste de modelos al variograma
modelos <- list("Exp", "Sph", "Gau", "Exc", "Mat", "Bes", "Pen", "Hol", "Wav","Pen+Wav")
lsaKs <- adjustMc(modelos,clasicVar,kappa=TRUE) # lsa = least square adjustment
mlKs <- adjustML(geoPreciTrans, modelos, geoPreciTrans$Err) #ml = ajuste x máxima verosimilitud
olsKs <- plotAdVarMc(clasicVar,lsaKs,TRUE,distVec=hSeq,data=preci.df)
wlsKs <- plotAdVarMc(clasicVar,lsaKs,FALSE,distVec=hSeq,data=preci.df)
cat("El mejor modelo por estadísticos en OLS: ",olsKs$resultados_OLS$Modelo[which.min(olsKs$resultados_OLS$BstStd[olsKs$resultados_OLS$BstStd != 0])],". Con un valor de: ", 
    olsKs$resultados_OLS$BstStd[which.min(olsKs$resultados_OLS$BstStd[olsKs$resultados_OLS$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en WLS: ",wlsKs$resultados_WLS$Modelo[which.min(wlsKs$resultados_WLS$BstStd[wlsKs$resultados_WLS$BstStd != 0])],". Con un valor de: ", 
    wlsKs$resultados_WLS$BstStd[which.min(wlsKs$resultados_WLS$BstStd[wlsKs$resultados_WLS$BstStd != 0])], "\n")
plotAdML(clasicVar, mlKs, method="ML")
plotAdML(clasicVar, mlKs, method="REML")
cat("El mejor modelo por estadísticos en ML: ",mlKs$resumen_ML$modelo[which.min(mlKs$resumen_ML$BstStd[mlKs$resumen_ML$BstStd != 0])],"con un valor de: ",
    mlKs$resumen_ML$BstStd[which.min(mlKs$resumen_ML$BstStd[mlKs$resumen_ML$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en REML: ",mlKs$resumen_REML$modelo[which.min(mlKs$resumen_REML$BstStd[mlKs$resumen_REML$BstStd != 0])],"con un valor de: ",
    mlKs$resumen_REML$BstStd[which.min(mlKs$resumen_REML$BstStd[mlKs$resumen_REML$BstStd != 0])], "\n")

# --- Interpolación Kriging Simple --- #
gridded(puntos) <- TRUE
modFinKs <- lsaKs$OLS$Gau
ks <- krige(formula=z ~ 1,locations=SpatialPointsDataFrame(coords=preci.trans[,c("x","y")],data=preci.trans,
                                                           proj4string=CRS("+init=epsg:3116")),
            newdata=puntos,model=modFinKs,beta=medValues$Ponderador[[5]])
predKs <- cbind(coordinates(ks), ks@data)
predDbKs <- db.create(predKs, autoname = FALSE)
predDbBtKs <- anam.y2z(predDbKs, names = "var1.pred", anam = mdb.herm)
coordinates(predKs) <- ~x1 + x2
gridded(predKs) <- TRUE
proj4string(predKs) <- CRS("+init=epsg:3116")
predKs$pred.Transform <- predDbBtKs@items$Raw.var1.pred
# --- Mapa de Predicción -- #
# Generar los objetos de gráfico spplot
mapPredKs<- spplot(predKs["pred.Transform"], cuts = 60, scales = list(draw = TRUE),
                    xlab = "Este (m)", ylab = "Norte (m)", main = "Predicción con Kriging Simple",
                    auto.key = FALSE)

mapVarKs <- spplot(predKs["var1.var"], cuts = 60, scales = list(draw = TRUE),
                    xlab = "Este (m)", ylab = "Norte (m)", main = "Varianza de la Predicción con Kriging Simple",
                    auto.key = FALSE)

# Mostrar ambos en una misma ventana
x11(width = 12, height = 6)
grid.arrange(mapPredKs, mapVarKs, ncol = 2)
#--------------------------------------------#
#             Kriging Ordinario              #
#--------------------------------------------#

# Verificación tendencia con datos transformados
summary(lm(z~x+y,data=preci.trans))

# Verificando anisotropia de los Errores
# estimateAnisotropy(SpatialPointsDataFrame(coords = preci.trans[, c("x", "y")], data = preci.trans,
#                         proj4string = CRS("+init=epsg:3116")),"z")
var4 <- variog4(geoPreciTrans, data = geoPreciTrans$data, max.dist=rango)
plot(var4, col = c("red", "blue", "green", "purple"), lwd = 2)

# Calculo  del Variograma Experimental para los datos tomados en campo
library(spatial)
library(sgeostat)

lags <- 15
preci.pts <- point(preci.trans, x="x", y ="y")
preci.pair <- pair(preci.pts, num.lags=lags,maxdist=rango)
preci.vars <- est.variograms(preci.pts,preci.pair,"z",trim=0.1)
ylim <- range(preci.vars[, c("robust", "med", "classic", "trimmed.mean")], na.rm = TRUE)
plot(preci.vars$bins,preci.vars$robust,lty=1, col=1,main = "Modelos Semivarianzas Empiricos",xlab="Distancia", ylab="Semivarianza", type="l",ylim=ylim)
lines(preci.vars$bins,preci.vars$med, col=2)
lines(preci.vars$bins,preci.vars$classic, col=3)
lines(preci.vars$bins,preci.vars$trimmed.mean, col=4)
legend("bottomright", c("Robusto", "Mediana", "Clásico","Media Recortada"), col=c(1,2,3,4), lty=c(1,1,1,1))

detach("package:spatial")
detach("package:sgeostat")

# Variograma Experimental con variogram
clasicVar <- variogram(z~1,~x+y,data = preci.trans, cutoff=rango)
cressieVar <- variogram(z~1,~x+y,data = preci.trans, cressie = TRUE,cutoff=rango)

# Variograma Experimental para Kriging Ordinario con est.variograms
koExp <- varExp(preci.vars,lags,"classic")
hSeq <- seq(0, max(koExp$dist, na.rm = TRUE), length.out = 15)
varExpKo <- variog(geodata=geoPreciTrans,coords=geoPreciTrans$coords,data=geoPreciTrans$data, 
                   trend="cte", uvec=hSeq,option="bin",max.dist=max(koExp$dist))

# Ploteo Ambos Variogramas
plot(koExp$dist, koExp$gamma, main = expression(paste("Comparación de ", gamma(h), " experimental")),
     xlab = expression("Distancia " ~ h), ylab = expression(gamma(h)), type="l",
     lwd=2, col="black",xlim=c(0,max(c(koExp$dist, varExpKo$u), na.rm = TRUE)),ylim=c(0,max(c(koExp$gamma, varExpKo$v), na.rm = TRUE)))
lines(varExpKo$u,varExpKo$v,col="red",type="l",lwd=2)
lines(clasicVar$dist, clasicVar$gamma, col="orange",type="l", lwd="2")
legend("bottomright", legend = c("Variograma con est.variograms", "Variogramas con variog", "Clásico con Variogram"),
       col = c("black", "red","orange"), lty = 1, lwd = 2)

# Ajuste de modelos al variograma
modelos <- list("Exp", "Sph", "Gau", "Exc", "Mat", "Bes", "Pen", "Hol", "Wav","Pen+Wav")
lsaKo <- adjustMc(modelos,clasicVar,kappa=TRUE) # lsa = least square adjustment
mlKo <- adjustML(geoPreciTrans, modelos, geoPreciTrans$data) #ml = ajuste x máxima verosimilitud
olsKo <- plotAdVarMc(clasicVar,lsaKo,TRUE,distVec=hSeq,data=preci.trans)
wlsKo <- plotAdVarMc(clasicVar,lsaKo,FALSE,distVec=hSeq,data=preci.trans)
cat("El mejor modelo por estadísticos en OLS: ",olsKo$resultados_OLS$Modelo[which.min(olsKo$resultados_OLS$BstStd[olsKo$resultados_OLS$BstStd != 0])],". Con un valor de: ", 
    olsKo$resultados_OLS$BstStd[which.min(olsKo$resultados_OLS$BstStd[olsKo$resultados_OLS$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en WLS: ",wlsKo$resultados_WLS$Modelo[which.min(wlsKo$resultados_WLS$BstStd[wlsKo$resultados_WLS$BstStd != 0])],". Con un valor de: ", 
    wlsKo$resultados_WLS$BstStd[which.min(wlsKo$resultados_WLS$BstStd[wlsKo$resultados_WLS$BstStd != 0])], "\n")
plotAdML(clasicVar, mlKo, method="ML")
plotAdML(clasicVar, mlKo, method="REML")
cat("El mejor modelo por estadísticos en ML: ",mlKo$resumen_ML$modelo[which.min(mlKo$resumen_ML$BstStd[mlKo$resumen_ML$BstStd != 0])],"con un valor de: ",
    mlKo$resumen_ML$BstStd[which.min(mlKo$resumen_ML$BstStd[mlKo$resumen_ML$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en REML: ",mlKo$resumen_REML$modelo[which.min(mlKo$resumen_REML$BstStd[mlKo$resumen_REML$BstStd != 0])],"con un valor de: ",
    mlKo$resumen_REML$BstStd[which.min(mlKo$resumen_REML$BstStd[mlKo$resumen_REML$BstStd != 0])], "\n")

# --- Interpolación Kriging Simple --- #
gridded(puntos) <- TRUE
modFinKo <- lsaKo$OLS$Gau
ko <- krige(formula=z ~ 1,locations=SpatialPointsDataFrame(coords=preci.trans[,c("x","y")],data=preci.trans,
                                                           proj4string=CRS("+init=epsg:3116")),
            newdata=puntos,model=modFinKo)
predKo <- cbind(coordinates(ko), ko@data)
predDbKo <- db.create(predKo, autoname = FALSE)
predDbBtKo <- anam.y2z(predDbKo, names = "var1.pred", anam = mdb.herm)
coordinates(predKo) <- ~x1 + x2
gridded(predKo) <- TRUE
proj4string(predKo) <- CRS("+init=epsg:3116")
predKo$pred.Transform <- predDbBtKo@items$Raw.var1.pred
# --- Mapa de Predicción -- #
# Generar los objetos de gráfico spplot
mapPredKo <- spplot(predKo["pred.Transform"], cuts = 60, scales = list(draw = TRUE),
                   xlab = "Este (m)", ylab = "Norte (m)", main = "Predicción con Kriging Ordinario",
                   auto.key = FALSE)

mapVarKo <- spplot(predKo["var1.var"], cuts = 60, scales = list(draw = TRUE),
                   xlab = "Este (m)", ylab = "Norte (m)", main = "Varianza de la Predicción con Kriging Ordinario",
                   auto.key = FALSE)

# Mostrar ambos en una misma ventana
x11(width = 12, height = 6)
grid.arrange(mapPredKo, mapVarKo, ncol = 2)

#--------------------------------------------#
#             Kriging Indicador              #
#--------------------------------------------#
# Valor crítico
preci.df$I <- ifelse(preci.df$z > 37, 1, 0)
geoPreci$I <- ifelse(geoPreci$data > 37, 1, 0)
# Verificación tendencia con datos transformados
summary(lm(I ~ x + y, data = preci.df))

# Verificación anisotropia
#estimateAnisotropy(SpatialPointsDataFrame(coords = preci.df[, c("x", "y")], data = preci.df,
#                                           proj4string = CRS("+init=epsg:3116")), "I")
var4 <- variog4(geoPreci, data = geoPreci$I, max.dist = rango)
plot(var4, col = c("red", "blue", "green", "purple"), lwd = 2)

# Calculo  del Variograma Experimental para Krigeado Indicador
library(spatial)
library(sgeostat)

lags <- 15
preci.pts <- point(preci.df, x="x", y ="y")
preci.pair <- pair(preci.pts, num.lags=lags,maxdist=rango)
preci.vars <- est.variograms(preci.pts,preci.pair,"I",trim=0.1)
ylim <- range(preci.vars[, c("robust", "classic", "trimmed.mean")], na.rm = TRUE)
plot(preci.vars$bins,preci.vars$robust,lty=1, col=1,main = "Modelos Semivarianzas Empiricos",xlab="Distancia", ylab="Semivarianza", type="l",ylim=ylim)
lines(preci.vars$bins,preci.vars$classic, col=2)
lines(preci.vars$bins,preci.vars$trimmed.mean, col=3)
legend("bottomright", c("Robusto", "Clásico","Media Recortada"), col=c(1,2,3), lty=c(1,1,1))

detach("package:spatial")
detach("package:sgeostat")

# Variograma Experimental con variogram
cressieVar <- variogram(I~1,~x+y,data = preci.df, cressie = TRUE,cutoff=rango)

# Variograma Experimental para Kriging Ordinario con est.variograms
kiExp <- varExp(preci.vars,lags,"robust")
hSeq <- seq(0, max(kiExp$dist, na.rm = TRUE), length.out = 15)
varExpKi <- variog(geodata=geoPreci,coords=geoPreciTrans$coords,data=geoPreci$I, 
                   trend="cte", uvec=hSeq,option="bin",max.dist=max(kiExp$dist), estimator.type = "modulus")

# Ploteo Ambos Variogramas
plot(kiExp$dist, kiExp$gamma, main = expression(paste("Comparación de ", gamma(h), " experimental")),
     xlab = expression("Distancia " ~ h), ylab = expression(gamma(h)), type="l",
     lwd=2, col="black",xlim=c(0,max(c(kiExp$dist, varExpKi$u), na.rm = TRUE)),ylim=c(0,max(c(kiExp$gamma, varExpKi$v), na.rm = TRUE)))
lines(varExpKi$u,varExpKi$v,col="red",type="l",lwd=2)
lines(cressieVar$dist, cressieVar$gamma, col="green",type="l", lwd="2")
legend("bottomright", legend = c("Variograma con est.variograms", "Variogramas con variog", "Robusto con Variogram"),
       col = c("black", "red","green"), lty = 1, lwd = 2)

# Ajuste de modelos al variograma
modelos <- list("Exp", "Sph", "Gau", "Exc", "Mat", "Bes", "Pen", "Hol", "Wav","Pen+Wav")
lsaKi <- adjustMc(modelos,cressieVar,kappa=TRUE) # lsa = least square adjustment
mlKi <- adjustML(geoPreci, modelos, geoPreci$I) #ml = ajuste x máxima verosimilitud
olsKi <- plotAdVarMc(cressieVar,lsaKi,TRUE,distVec=hSeq,data=preci.df)
wlsKi <- plotAdVarMc(cressieVar,lsaKi,FALSE,distVec=hSeq,data=preci.df)
cat("El mejor modelo por estadísticos en OLS: ",olsKi$resultados_OLS$Modelo[which.min(olsKi$resultados_OLS$BstStd[olsKi$resultados_OLS$BstStd != 0])],". Con un valor de: ", 
    olsKi$resultados_OLS$BstStd[which.min(olsKi$resultados_OLS$BstStd[olsKi$resultados_OLS$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en WLS: ",wlsKi$resultados_WLS$Modelo[which.min(wlsKi$resultados_WLS$BstStd[wlsKi$resultados_WLS$BstStd != 0])],". Con un valor de: ", 
    wlsKi$resultados_WLS$BstStd[which.min(wlsKi$resultados_WLS$BstStd[wlsKi$resultados_WLS$BstStd != 0])], "\n")
plotAdML(cressieVar, mlKi, method="ML")
plotAdML(cressieVar, mlKi, method="REML")
cat("El mejor modelo por estadísticos en ML: ",mlKi$resumen_ML$modelo[which.min(mlKi$resumen_ML$BstStd[mlKi$resumen_ML$BstStd != 0])],"con un valor de: ",
    mlKi$resumen_ML$BstStd[which.min(mlKi$resumen_ML$BstStd[mlKi$resumen_ML$BstStd != 0])], "\n")
cat("El mejor modelo por estadísticos en REML: ",mlKi$resumen_REML$modelo[which.min(mlKi$resumen_REML$BstStd[mlKi$resumen_REML$BstStd != 0])],"con un valor de: ",
    mlKi$resumen_REML$BstStd[which.min(mlKi$resumen_REML$BstStd[mlKi$resumen_REML$BstStd != 0])], "\n")

# --- Interpolación Kriging Indicador --- #

gridded(puntos) <- TRUE
modFinKi <- lsaKi$OLS$Hol
ki <- krige(I(preci.df$z>37)~1, locations = SpatialPointsDataFrame( coords = preci.df[, c("x", "y")], data = preci.df,
    proj4string = CRS("+init=epsg:3116")),newdata = puntos, model = modFinKi)
predKi <- cbind(coordinates(ki), ki@data)
coordinates(predKi) <- ~ x1 + x2
gridded(predKi) <- TRUE
proj4string(predKi) <- CRS("+init=epsg:3116")
# --- Mapa de Predicción -- #
# Generar los objetos de gráfico spplot
mapPredKi <- spplot(predKi["var1.pred"],
                    cuts = 60, scales = list(draw = TRUE),
                    xlab = "Este (m)", ylab = "Norte (m)", main = "Predicción con Kriging Indicador",
                    auto.key = FALSE
)

mapVarKi <- spplot(predKi["var1.var"],
                   cuts = 60, scales = list(draw = TRUE),
                   xlab = "Este (m)", ylab = "Norte (m)", main = "Varianza de la Predicción con Kriging Indicador",
                   auto.key = FALSE
)

# Mostrar ambos en una misma ventana
x11(width = 12, height = 6)
grid.arrange(mapPredKi, mapVarKi, ncol = 2)

#--------------------------------------------#
#             Validación Cruzada             #
#--------------------------------------------#

ksCvZ <- krige.cv(
  formula = z ~ 1, locations = SpatialPointsDataFrame(coords = preci.trans[, c("x", "y")], data = preci.trans, proj4string = CRS("+init=epsg:3116")),
  model = modFinKs, beta = medValues$Ponderador[[5]], nmax = 15
)
koCvZ <- krige.cv(
  formula = z ~ 1, locations = SpatialPointsDataFrame(coords = preci.trans[, c("x", "y")], data = preci.trans, proj4string = CRS("+init=epsg:3116")),
  model = modFinKo, nmax = 15
)

KiCvI <- krige.cv(
  formula = I ~ 1, locations = SpatialPointsDataFrame(coords = preci.df[, c("x", "y")], data = preci.df, proj4string = CRS("+init=epsg:3116")),
  model = modFinKi, nmax=15
)

# Concordantes para KI
concor1 <- nrow(KiCvI[KiCvI$var1.pred>=0.5 & KiCvI$var1.observed>0])
concor0 <- nrow(KiCvI[KiCvI$var1.pred<0.5 & KiCvI$var1.observed<1])
conc <- concor1 + concor0
discor10 <- nrow(KiCvI[KiCvI$var1.pred>=0.5 & KiCvI$var1.observed<1])
discor01 <- nrow(KiCvI[KiCvI$var1.pred<0.5 & KiCvI$var1.observed>0])
discor <- discor10 + discor01
empates <- sum(ifelse(KiCvI$var.pred==0.5,1,0))
gamma <- (conc - discor)/(conc + discor + empates)

resCvZ <- rbind(criteria.cv(ksCvZ), criteria.cv(koCvZ),criteria.cv(KiCvI))
rownames(resCvZ) <- c("Kriging Simple", "Kriging Ordinario", "Kriging Indicador")
resCvZ

##############################################
#                Diseño Redes                #
##############################################

cundi.sf <- st_as_sf(cundi.sp)
borde <- st_boundary(st_union(cundi.sf))
cords <- st_coordinates(borde)
p1 <- Polygon(cords[,1:2])
ps1 <- Polygons(list(p1), "1")
poly <- SpatialPolygons(list(ps1))

netPreci.ks1<- network.design(z~1, modFinKs, npoints=100, boundary=poly, nmax=15, type="random", beta=medValues$Ponderador[[5]])
netPreci.ks2<- network.design(z~1, modFinKs, npoints=400, boundary=poly, nmax=15, type="random", beta=medValues$Ponderador[[5]])
netPreci.ks3<- network.design(z~1, modFinKs, npoints=900, boundary=poly, nmax=15, type="random", beta=medValues$Ponderador[[5]])
netPreci.ks4<- network.design(z~1, modFinKs, npoints=1600, boundary=poly, nmax=15, type="random", beta=medValues$Ponderador[[5]])
netPreci.ks5<- network.design(z~1, modFinKs, npoints=2500, boundary=poly, nmax=15, type="random", beta=medValues$Ponderador[[5]])

netPreci.ko1 <- network.design(z~1, modFinKo, npoints=100, boundary=poly, nmax=15, type="random")
netPreci.ko2 <- network.design(z~1, modFinKo, npoints=400, boundary=poly, nmax=15, type="random")
netPreci.ko3 <- network.design(z~1, modFinKo, npoints=900, boundary=poly, nmax=15, type="random")
netPreci.ko4 <- network.design(z~1, modFinKo, npoints=1600, boundary=poly, nmax=15, type="random")
netPreci.ko5 <- network.design(z~1, modFinKo, npoints=2500, boundary=poly, nmax=15, type="random")

