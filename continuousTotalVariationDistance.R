continuousTotalVariationDistance <- function (x, y, w_x = NULL, w_y = NULL,
                                              lower = -Inf, upper = Inf,
                                              tol = 1e-5,
                                              ...) {
  #' continuousTotalVariationDistance
  #'
  #' Computes a 1D Total Variation distances (TVd) estimate between two classes for continuos data. Based on density estimation and numerical integration.
  #'
  #' @param x A 1D array.
  #' @param y A 1D array.
  #' @param w_x Weights for x, will be normalised to one.
  #' @param w_y Weights for x, will be normalised to one.
  #' @param lower Lower limit for the numerical integration. 
  #' @param upper Upper limit for the numerical integration.
  #' @param tol Tolerance for the numerical integration.
  #' @param ... Parameters for the density estimation. 
  #' 
  #' @return An estimate of the total variation distance.  
  #' 
  #' @examples
  #' data.example <- data.frame(V1 = rnorm(150), V2 = rnorm(150, 0.5, 1),
  #'  Treatment = c(rep(0, 100), rep(1, 50)),
  #'  Hospitalisation = c(rep(0,70), rep(1, 30), rep(1, 20), rep(0, 30)),
  #'  Death = (runif(150) < 0.1) * 1)
  #' 
  #' TVd <- continuousTotalVariationDistance(data.example$V1[data.example$Treatment == 0],
  #'  data.example$V1[data.example$Treatment == 1], lower = -5, upper = 5)
  #'  
  #' @references H Inouzhe, I Barrio, MX Rodríguez-Álvarez, P Gordaliza, I Bengoechea and JM Quintana. (2023) A study on group fairness in healthcare outcomes for nursing home residents during the COVID-19 pandemic in the Basque Country.
  #'
  #' @export
  #'
  
  if(is.null(w_x)){
    w_x <- rep(1, length(x))
  } else {
    if(sum(w_x) != 1){
      warning("w_x does not sum to 1. We will normalise it!")
      w_x <- w_x / sum(w_x)
    }
  }
  if(is.null(w_y)){
    w_y <- rep(1, length(y))
  } else {
    if(sum(w_y) != 1){
      warning("w_y does not sum to 1. We will normalise it!")
      w_y <- w_y / sum(w_y)
    }
  }
  fx_0 <- tryCatch(ks::kde (as.numeric(x), w = w_x, ...), error = function(e){
    bw <- density(as.numeric(x))$bw
    return(ks::kde (as.numeric(x), h = bw, w = w_x, ...))
  })
  fy_0 <- tryCatch(ks::kde (as.numeric(y), w = w_y, ...), error = function(e){
    bw <- density(as.numeric(y))$bw
    return(ks::kde (as.numeric(y), h = bw, w = w_y, ...))
  })
  fx <- function (t) predict(fx_0, x = t) 
  fy <- function (t) predict(fy_0, x = t) 
  
  g <- function(z) abs(fx(z) - fy(z))
  cubature::hcubature(g, lower, upper, tol = tol)$integral / 2
}
docstring::docstring(continuousTotalVariationDistance)