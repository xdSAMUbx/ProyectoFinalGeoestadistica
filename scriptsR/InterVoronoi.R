library(deldir)
library(RGeostats)
library(geoR)
library(gstat)
library(spatial)
library(geospt)
library(geosptdb)
library(sgeostat)
library(sp)
library(sf)
library(intamap)
library(FNN)
library(readxl)
library(scatterplot3d)
library(tidyr)
library(lattice)
library(dplyr)
library(pbapply)

####################
# Cargando la Data #
####################

rm(list = ls())
setwd("D:\\Universidad\\7 Semestre\\Geoestadistica\\Ejercicios\\Datos de Area\\EjercicioSAR\\proyectoFinal")
path <- "../proyectoFinal/PluvidbFinal.csv"
pathCundi <- "../proyectoFinal/munCundiShp/Municipios_Cundinamarca.shp"
# pathCundi <- "D:/Universidad/7 semestre/Geoestadistica/ProyectoFinal/Pluviometria/munCundiShp/Municipios_Cundinamarca.shp"
source("../proyectoFinal/polVor.R")
save.image(file="../proyectoFinal/.Rdata")
load("../proyectoFinal/.Rdata")
data <- read.csv(path)

cundi <- st_read(pathCundi)
cundi <- st_transform(cundi, crs = 3116)

df <- dplyr::select(data, Nombre, Latitud, Longitud, Municipio, starts_with("X"))
for (col in names(df)[5:ncol(df)]) {
  df[[col]] <- as.numeric(gsub(",", ".", df[[col]]))
}

preci <- dplyr::select(df, Nombre, Latitud,Longitud,Municipio,"X2020.01")
preci[preci == 0] <- NA
preci <- na.omit(preci)
preci<- st_as_sf(preci, coords = c("Longitud", "Latitud"), crs = 4326)
preci <- st_transform(preci, crs = 3116)
preci <- preci %>% mutate(Municipio = toupper(Municipio))
preci <- preci %>% rename("z" = "X2020.01")

x11()
plot(st_geometry(cundi), col = "lightgray", border = "black", main = "Mapa de Cundinamarca")
plot(st_geometry(preci), col = "red", pch = 20, cex = 1.2, add = TRUE)

vor <- polVor(preci,"Nombre",knn=20,TRUE,TRUE,cundi$geometry,TRUE)
plot(vor$geometry, border = "red")
plot(st_geometry(preci), add = TRUE, pch = 19, col = "blue", cex=0.3)

cundiSp <- as(cundi,"Spatial")
pts <- spsample(cundiSp,n=30000,type="regular")
coordsCundi <- coordinates(pts)
colnames(coordsCundi) <- c("x","y")
coordsPred <- as.data.frame(coordsCundi)
coordsSf <- st_as_sf(coordsPred, coords = c("x", "y"), crs = 3116)

library(parallel)

interpolar_voronoi_masivo <- function(
    puntos_nuevos,
    estaciones_sf,
    variable = "valor",
    n.cores = detectCores() - 1
) {
  stopifnot(inherits(puntos_nuevos, "sf"),
            inherits(estaciones_sf, "sf"),
            variable %in% names(estaciones_sf))
  
  # 1) Aseguro columna 'id' en estaciones
  if (!"id" %in% names(estaciones_sf)) {
    estaciones_sf$id <- paste0("est_", seq_len(nrow(estaciones_sf)))
  }
  
  # 2) Genero UN ÚNICO voronoi de las estaciones
  message(" ▸ Generando Voronoi de las estaciones (una sola vez)…")
  vor_geom_est <- st_voronoi(st_union(estaciones_sf))
  vor_sfc_est  <- st_collection_extract(
    st_sfc(vor_geom_est, crs = st_crs(estaciones_sf)),
    "POLYGON"
  )
  voronoi_est <- st_sf(id = estaciones_sf$id, geometry = vor_sfc_est)
  
  # 3) Preparo cluster
  cl <- makeCluster(n.cores)
  clusterEvalQ(cl, { library(sf) })
  clusterExport(
    cl,
    varlist = c("puntos_nuevos", "estaciones_sf", "voronoi_est", "variable"),
    envir   = environment()
  )
  
  message(sprintf(
    " ▸ Interpolando %d puntos usando %d núcleos…",
    nrow(puntos_nuevos), n.cores
  ))
  
  # 4) Loop paralelo, uno por punto
  resultados <- parLapply(cl, seq_len(nrow(puntos_nuevos)), function(i) {
    p_nuevo <- puntos_nuevos[i, , drop = FALSE]
    
    # 4.1) Voronoi temporal para estaciones + este punto
    geoms       <- c(st_geometry(estaciones_sf), st_geometry(p_nuevo))
    combo_sf    <- st_sf(geometry = st_sfc(geoms, crs = st_crs(estaciones_sf)))
    vor_tmp     <- st_voronoi(st_union(combo_sf))
    vor_sfc_tmp <- st_collection_extract(
      st_sfc(vor_tmp, crs = st_crs(estaciones_sf)),
      "POLYGON"
    )
    vor_tmp_sf  <- st_sf(geometry = vor_sfc_tmp)
    
    # 4.2) Me quedo con el último polígono (el del punto nuevo)
    pol_nuevo <- vor_tmp_sf[nrow(vor_tmp_sf), , drop = FALSE]
    
    # 4.3) Intersección con los Voronoi de estaciones
    inters <- suppressWarnings(st_intersection(pol_nuevo, voronoi_est))
    if (nrow(inters) == 0) return(NA_real_)
    
    # 4.4) Calculo áreas e pesos normalizados
    inters$area <- st_area(inters)
    total_area <- sum(inters$area)
    inters$peso <- inters$area / total_area
    
    # 4.5) Extraigo valores de la estación original y sumo ponderado
    vals <- estaciones_sf[[variable]][ match(inters$id, estaciones_sf$id) ]
    sum(vals * inters$peso, na.rm = TRUE)
  })
  
  stopCluster(cl)
  message(" ▸ Interpolación completada.")
  
  unlist(resultados)
}
prueba <- interpolar_voronoi_masivo(coordsSf,preci,vor$data,"z")
# Extraer coordenadas
coords <- st_coordinates(df)
df$X <- coords[,1]
df$Y <- coords[,2]
df <- as.data.frame(df)  # Salimos de sf para usar con geoR
