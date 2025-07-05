varExp <- function(objeto, lags, tipo) {
  stopifnot(inherits(objeto, "variogram"))
  stopifnot(is.character(tipo) && length(tipo) == 1)
  
  dir.hor <- rep(0, lags)
  dir.ver <- rep(0, lags)
  id <- rep("var1", lags)
  
  y <- data.frame(
    np = objeto$n[1:lags],
    dist = objeto$bins[1:lags],
    gamma = objeto[[tipo]][1:lags],
    dir.hor = dir.hor,
    dir.ver = dir.ver,
    id = id
  )
  
  class(y) <- c("variogram", "gstatVariogram", "data.frame")
  return(y)
}