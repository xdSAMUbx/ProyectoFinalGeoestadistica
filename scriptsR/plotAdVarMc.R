plotAdVarMc <- function(varExp, adj=NULL, mco=TRUE, distVec, data) {
  stopifnot(inherits(varExp, "gstatVariogram") | inherits(varExp, "variogram"))
  stopifnot(!is.null(adj))
  
  # Data frame para almacenar resultados
  resultados_OLS <- data.frame(
    Modelo = character(),
    MPE = numeric(),
    RMSPE = numeric(),
    R2 = numeric(),
    BstStd = numeric(),
    stringsAsFactors = FALSE
  )
  
  resultados_WLS <- data.frame(
    Modelo = character(),
    MPE = numeric(),
    RMSPE = numeric(),
    R2 = numeric(),
    BstStd = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Ajuste por Mínimos Cuadrados Ordinarios (MCO)
  if (mco == TRUE) {
    plot(varExp$dist, varExp$gamma, pch = 20, col = "black", type = "p", 
         xlab = "Distancia", ylab = "Semivarianza", 
         main = "Ajuste por Mínimos Cuadrados Ordinarios (MCO)")
    
    colores <- rainbow(length(adj$OLS))
    i <- 1
    
    for (nombre in names(adj$OLS)) {
      modelo <- adj$OLS[[nombre]]
      if (inherits(modelo, "try-error") || is.null(modelo)) next
      lines(variogramLine(modelo, maxdist = max(varExp$dist)), col = colores[i], lwd = 2)
      linea <- variogramLine(modelo, dist_vector = distVec,col = colores[i], lwd = 2)
      # Cálculos de estadísticos
      mpe <- sum(varExp$gamma - linea$gamma) / nrow(data)
      rmspe <- sqrt(sum((varExp$gamma - linea$gamma)^2) / nrow(data))
      meanz <- mean(varExp$gamma)
      r2 <- 1 - sum((linea$gamma - varExp$gamma)^2) / sum((varExp$gamma - meanz)^2)
      BstStd <- sqrt(mpe^2 + rmspe^2 + (r2 - 1)^2)
      
      # Almacenar resultados en el data frame
      resultados_OLS <- rbind(resultados_OLS, data.frame(
        Modelo = nombre, MPE = mpe, RMSPE = rmspe, R2 = r2, BstStd = BstStd
      ))
      
      i <- i + 1
    }
    
    legend("bottomright", legend = names(adj$OLS), col = colores, lty = 1, lwd = 2, cex = 1)
    
  } else {  # Ajuste por Mínimos Cuadrados Ponderados (MCP)
    plot(varExp$dist, varExp$gamma, pch = 20, col = "black", type = "p", 
         xlab = "Distancia", ylab = "Semivarianza", 
         main = "Ajuste por Mínimos Cuadrados Ponderados (MCP)")
    
    colores <- rainbow(length(adj$WLS))
    i <- 1
    
    for (nombre in names(adj$WLS)) {
      modelo <- adj$WLS[[nombre]]
      if (inherits(modelo, "try-error") || is.null(modelo)) next
      lines(variogramLine(modelo, maxdist = max(varExp$dist)), col = colores[i], lwd = 2)
      linea <- variogramLine(modelo, dist_vector = distVec,col = colores[i], lwd = 2)
      
      # Cálculos de estadísticos
      mpe <- sum(varExp$gamma - linea$gamma) / nrow(data)
      rmspe <- sqrt(sum((varExp$gamma - linea$gamma)^2) / nrow(data))
      meanz <- mean(varExp$gamma)
      r2 <- 1 - sum((linea$gamma - varExp$gamma)^2) / sum((varExp$gamma - meanz)^2)
      BstStd <- sqrt(mpe^2 + rmspe^2 + (r2 - 1)^2)
      
      # Almacenar resultados en el data frame
      resultados_WLS <- rbind(resultados_WLS, data.frame(
        Modelo = nombre, MPE = mpe, RMSPE = rmspe, R2 = r2, BstStd = BstStd
      ))
      
      i <- i + 1
    }
    
    legend("bottomright", legend = names(adj$WLS), col = colores, lty = 1, lwd = 2, cex = 1)
  }
  
  # Devolver los data frames con los resultados
  return(list(resultados_OLS = resultados_OLS, resultados_WLS = resultados_WLS))
}
