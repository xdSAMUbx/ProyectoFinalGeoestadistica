polVor <- function(objeto, id=NULL, knn = 10, n.neig = FALSE, area = FALSE, mask = NULL, paralel = FALSE) {
  stopifnot(inherits(objeto, "sf"))
  
  coords <- sf::st_coordinates(objeto)
  vor <- deldir::tile.list(deldir::deldir(coords[,1], coords[,2]))
  
  if (is.null(id)) {
    objeto$id_temp__ <- seq_len(nrow(objeto))
    id <- "id_temp__"
  }
  
  # Crear polígonos (secuencial o paralelo)
  if (paralel) {
    # Usar un núcleo menos que el total disponible
    n.cores <- max(1, parallel::detectCores() - 1)
    cl <- parallel::makeCluster(n.cores)
    
    # Exportar variables y cargar librerías en los workers
    parallel::clusterExport(cl, varlist = c("vor"), envir = environment())
    parallel::clusterEvalQ(cl, {
      library(sp)
    })
    
    # Usar pbapply con barra de progreso
    vor_polygons <- pbapply::pblapply(vor, function(tile) {
      sp::Polygons(list(sp::Polygon(cbind(tile$x, tile$y))), ID = as.character(tile$ptNum))
    }, cl = cl)
    
    parallel::stopCluster(cl)
    message(paste("Paralelización finalizada con", n.cores, "núcleos."))
    
  } else {
    # Secuencial con barra de progreso
    vor_polygons <- pbapply::pblapply(vor, function(tile) {
      sp::Polygons(list(sp::Polygon(cbind(tile$x, tile$y))), ID = as.character(tile$ptNum))
    })
  }
  
  # Crear SpatialPolygons y luego sf
  pol <- sf::st_as_sf(sp::SpatialPolygons(vor_polygons,
                                          proj4string = sp::CRS(sf::st_crs(objeto)$wkt)))
  pol$id <- objeto[[id]]
  
  # Aplicar máscara si se proporciona
  if (!is.null(mask)) {
    mask <- sf::st_transform(mask, sf::st_crs(pol))
    mask <- sf::st_union(mask)
    pol <- sf::st_intersection(pol, mask)
  }
  
  # Calcular pesos
  if (n.neig & !area) {
    d <- rowSums(FNN::get.knn(coords, knn)$nn.dist)
    pol$Wi <- d / sum(d)
    
  } else if (area & !n.neig) {
    a <- sf::st_area(pol)
    pol$Wi <- a / sum(a)
    
  } else if (area & n.neig) {
    d <- rowSums(FNN::get.knn(coords, knn)$nn.dist)
    a <- sf::st_area(pol)
    pol$WiNeigh <- d / sum(d)
    pol$WiArea <- a / sum(a)
    
  } else {
    message("No se calculó ningún peso.")
  }
  
  return(list(data = pol, geometry = sf::st_geometry(pol)))
}