meanDeclu <- function(objeto, dfs, knn = 10, nombres = NULL) {
  stopifnot(inherits(objeto, "matrix"))
  
  d <- rowSums(FNN::get.knn(objeto, knn)$nn.dist)
  Wi <- d / sum(d)
  
  pond_vals <- numeric(length(dfs))
  dfs_mod <- vector("list", length(dfs))
  
  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    colname <- intersect(c("z", "data"), names(df))[1]
    if (is.na(colname)) stop("No se encontró columna 'z' ni 'data'")
    
    z <- if (colname == "data") df$data else as.numeric(unlist(df[[colname]]))
    
    pond <- as.numeric(Wi %*% z)
    pond_vals[i] <- pond
    
    df$Wi <- Wi
    df$h <- d
    df$Err <- z - pond
    attr(df, "pond") <- pond
    
    dfs_mod[[i]] <- df
  }
  
  # Si no se pasan nombres, usar índices
  if (is.null(nombres)) {
    nombres <- paste0("obj", seq_along(dfs))
  }
  
  tabla_pond <- data.frame(
    Objeto = nombres,
    Ponderador = pond_vals
  )
  
  return(list(dfs = dfs_mod, tabla = tabla_pond))
}