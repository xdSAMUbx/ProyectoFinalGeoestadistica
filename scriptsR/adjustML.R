adjustML <- function(geodata, lista, data) {
  # Modelos disponibles y sus nombres en likfit
  mods <- list("Exp" = "exponential", "Sph" = "spherical", "Gau" = "gaussian",
               "Mat" = "matern", "Wav" = "wave", "Pow" = "power",
               "Exc" = "powered.exponential", "Nug" = "pure.nugget")
  
  if (is.list(lista)) lista <- unlist(lista)
  ML <- list(); REML <- list()
  
  # Modelos que usan kappa
  modelos_kappa <- c("Mat", "Exc")
  
  for (mod in lista) {
    if (mod %in% names(mods)) {
      usa_kappa <- mod %in% modelos_kappa
      
      cat("Ajustando por ML", mod, "con cov.model:", mods[[mod]], "\n")
      ML[[mod]] <- try(likfit(geodata, coords = geodata$coords, data = data, ini.cov.pars = c(10, 20),
                              cov.model = mods[[mod]], fix.nugget = FALSE, fix.kappa = !usa_kappa, 
                              trend = "cte",lik.method = "ML"),silent = TRUE)
      
      cat("Ajustando por REML", mod, "con cov.model:", mods[[mod]], "\n")
      REML[[mod]] <- try(likfit(geodata, coords = geodata$coords, data = data, ini.cov.pars = c(10, 20),
                                cov.model = mods[[mod]], fix.nugget = FALSE, fix.kappa = !usa_kappa,
                                trend = "cte", lik.method = "REML"), silent = TRUE)
    } else {
      cat(strrep("-", getOption("width")), "\n")
      cat("El modelo", mod, "no se puede ajustar por ML o REML\n")
      cat(strrep("-", getOption("width")), "\n")
    }
  }
  
  resumen_variogram <- function(lista_modelos) {
    modelos <- names(lista_modelos)
    df <- data.frame(
      modelo = modelos,
      meseta = NA_real_,
      rango  = NA_real_,
      kappa  = NA_real_,
      AIC    = NA_real_,
      BIC    = NA_real_,
      BstStd = NA_real_,
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(modelos)) {
      fit <- lista_modelos[[i]]
      if (inherits(fit, "try-error") || is.null(fit)) {
        # Si falló el ajuste, dejamos ceros
        df[i, c("meseta", "rango", "kappa", "AIC", "BIC", "BstStd")] <- 0
      } else {
        # fit$cov.pars == c(partial sill, range)
        df$meseta[i] <- fit$cov.pars[1]
        df$rango[i]  <- fit$cov.pars[2]
        
        # fit$kappa existe solo en modelos con kappa; si no, lo ponemos a 0
        df$kappa[i]  <- if (!is.null(fit$kappa)) fit$kappa else 0
        
        # Extraemos AIC y BIC
        df$AIC[i]    <- fit$AIC
        df$BIC[i]    <- fit$BIC
        
        # Cálculo del BstStd como la raíz de la suma de los cuadrados de AIC y BIC
        df$BstStd[i] <- sqrt(fit$AIC^2 + fit$BIC^2)
      }
    }
    return(df)
  }
  
  resumen_ML <- resumen_variogram(ML)
  resumen_REML <- resumen_variogram(REML)
  
  return(list(
    ML = ML,
    REML = REML,
    resumen_ML = resumen_ML,
    resumen_REML = resumen_REML
  ))
}