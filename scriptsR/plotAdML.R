plotAdML <- function(varExp, adj, method="ML") {
  stopifnot(inherits(varExp, "variogram") | inherits(varExp, "gstatVariogram"))
  stopifnot(method %in% c("ML", "REML"))
  
  # Selección del modelo según el método
  modelos <- if (method == "ML") adj$ML else adj$REML
  titulo <- if (method == "ML") "Ajuste por Máxima Verosimilitud (ML)" else "Ajuste por Máxima Verosimilitud Restringida (REML)"
  
  # Cálculo de los límites de la gráfica
  xlim <- c(0, max(varExp$dist, na.rm = TRUE))  # Rango de X según el variograma experimental
  ymax_exp <- max(varExp$gamma, na.rm = TRUE)   # Máximo valor de semivarianza experimental
  
  # Máximo de nugget + sigmasq de los modelos ajustados
  ymax_modelos <- max(sapply(modelos, function(m) {
    if (inherits(m, "try-error") || is.null(m)) return(NA)
    nug <- try(m[["nugget"]], silent = TRUE)
    sig <- try(m[["sigmasq"]], silent = TRUE)
    if (inherits(nug, "try-error") || inherits(sig, "try-error")) return(NA)
    nug + sig
  }), na.rm = TRUE)
  
  ylim <- c(0, max(ymax_exp, ymax_modelos, na.rm = TRUE))  # Rango de Y para la gráfica
  
  # Graficar variograma experimental con los límites calculados
  plot(varExp$dist, varExp$gamma, pch = 20, col = "black", type = "p",
       xlab = "Distancia", ylab = "Semivarianza", main = titulo, xlim = xlim, ylim = ylim)
  
  # Lista con modelos equivalentes
  mods <- list("Exp" = "exponential", "Sph" = "spherical", "Gau" = "gaussian",
               "Mat" = "matern", "Wav" = "wave", "Pow" = "power",
               "Exc" = "powered.exponential", "Nug" = "pure.nugget")
  
  # Definir colores para los modelos ajustados
  colores <- RColorBrewer::brewer.pal(length(modelos), "Set1")  # Usamos una paleta de colores fija
  
  # Añadir los modelos ajustados al gráfico
  i <- 1
  for (nombre in names(modelos)) {
    modelo <- modelos[[nombre]]
    
    if (inherits(modelo, "try-error") || is.null(modelo) || modelo$cov.pars[2]<=0) next  # Si el modelo falla, saltamos
    if (!is.null(modelo$cov.pars) && !is.null(modelo$cov.model)) {
      # Graficar la línea del modelo ajustado
      line_modelo <- variogramLine(vgm(modelo$cov.pars[1], names(mods)[which(mods == modelo$cov.model)], 
                                       modelo$cov.pars[2]), maxdist = max(varExp$dist))
      lines(line_modelo, col = colores[i], lwd = 2)
    }
    
    i <- i + 1
  }
  # Añadir la leyenda
  legend("bottomright", legend = names(modelos), col = colores, lty = 1, lwd = 2, cex = 0.9)
}