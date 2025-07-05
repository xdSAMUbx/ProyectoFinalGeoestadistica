# lista: lista de modelos que se desea ajustar
# ani: si hay presencia de variogramas anidados (máximo 2)
# variExp: Variograma Experimental
# kappa: si se tiene algún modelo dentro de la lista que necesite kappa

adjustMc <- function(lista, variExp, kappa = FALSE) {
  stopifnot(inherits(lista, c("list", "character")))
  stopifnot(inherits(variExp,"gstatVariogram") | (inherits(variExp,"variogram")))
  if (is.list(lista)) lista <- unlist(lista)
  modsKappa <- c("Mat", "Exc")
  OLS <- list(); WLS <- list()
  for (mod in lista) {
    if (grepl("^[A-Za-z0-9]+\\+[A-Za-z0-9]+$", mod)) {
      ani <- strsplit(mod, "\\+")[[1]]
      k1 <- ani[[1]] %in% modsKappa
      k2 <- ani[[2]] %in% modsKappa
      if (k2 || k1) {
        cat(strrep("-", getOption("width")), "\n")
        cat("Ajustando modelo Anidado (k)", mod, "por OLS\n")
        OLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, ani[[1]], NA, kappa = 0, add.to = vgm(NA, ani[[2]], NA, 0, kappa = 0)), fit.method = 6, fit.kappa = k1 | k2), silent = TRUE)
        cat("Ajustando modelok Anidado (k)", mod, "por WLS\n")
        WLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, ani[[1]], NA, kappa = 0, add.to = vgm(NA, ani[[2]], NA, 0, kappa = 0)), fit.method = 7, fit.kappa = k1 | k2), silent = TRUE)
      } else {
        cat(strrep("-", getOption("width")), "\n")
        cat("Ajustando modelo Anidado", mod, "por OLS\n")
        OLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, ani[[1]], NA, add.to = vgm(NA, ani[[2]], NA, 0)), fit.method = 6), silent = TRUE)
        cat("Ajustando modelo Anidado", mod, "por WLS\n")
        WLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, ani[[1]], NA, add.to = vgm(NA, ani[[2]], NA, 0)), fit.method = 7), silent = TRUE)
      }
    } else {
      if (mod %in% modsKappa) {
        cat(strrep("-", getOption("width")), "\n")
        cat("Ajustando modelo (k)", mod, "por OLS\n")
        OLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, mod, NA, 0, kappa = 0), fit.method = 6, fit.kappa = TRUE), silent = TRUE)
        cat("Ajustando modelo (k)", mod, "por WLS\n")
        WLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, mod, NA, 0, kappa = 0), fit.method = 7, fit.kappa = TRUE), silent = TRUE)
      } else {
        cat(strrep("-", getOption("width")), "\n")
        cat("Ajustando modelo", mod, "por OLS\n")
        OLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, mod, NA, 0), fit.method = 6), silent = TRUE)
        cat("Ajustando modelo", mod, "por WLS\n")
        WLS[[mod]] <- try(fit.variogram(variExp, vgm(NA, mod, NA, 0), fit.method = 7), silent = TRUE)
      }
    }
  }
  
  resumen_variogram <- function(lista_modelos) {
    nombres <- names(lista_modelos)
    df <- data.frame(
      modelo = nombres,
      anidado = grepl("\\+", nombres),
      nugget = numeric(length(nombres)),
      meseta1 = numeric(length(nombres)),
      rango1 = numeric(length(nombres)),
      meseta2 = numeric(length(nombres)),
      rango2 = numeric(length(nombres)),
      kappa1 = numeric(length(nombres)),
      kappa2 = numeric(length(nombres)),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(nombres)) {
      resultado <- lista_modelos[[i]]
      if (inherits(resultado, "try-error") || is.null(resultado)) {
        df[i, 3:9] <- 0
      } else {
        if (nrow(resultado) == 2) {
          df$nugget[i]  <- resultado$psill[1]
          df$meseta1[i] <- resultado$psill[2]
          df$rango1[i]  <- resultado$range[2]
          df$kappa1[i]  <- tryCatch(resultado$kappa[2], error = function(e) 0)
        } else if (nrow(resultado) >= 3) {
          df$nugget[i]  <- resultado$psill[1]
          df$meseta1[i] <- resultado$psill[2]
          df$rango1[i]  <- resultado$range[2]
          df$meseta2[i] <- resultado$psill[3]
          df$rango2[i]  <- resultado$range[3]
          df$kappa1[i]  <- tryCatch(resultado$kappa[2], error = function(e) 0)
          df$kappa2[i]  <- tryCatch(resultado$kappa[3], error = function(e) 0)
        }
      }
    }
    return(df)
  }
  
  dfOls <- resumen_variogram(OLS)
  dfWls <- resumen_variogram(WLS)
  
  return(list(OLS = OLS, WLS = WLS, dfOls = dfOls, dfWls = dfWls))
}
