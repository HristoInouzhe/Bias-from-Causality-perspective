weightedMeansDifference <- function(X, Y, weights.X, weights.Y, confidence.level = 0.05){
  #' weightedMeansDifference
  #'
  #' Computes a difference in weighted means and its confidence interval for one dimensional data.
  #'
  #' @param X A 1d array.
  #' @param Y A 1d array.
  #' @param weights.X Weights corresponding to X. No need to sum one.
  #' @param weights.Y Weights corresponding to X. No need to sum one.
  #' @param confidence.level Level for the confidence interval.
  #' 
  #' @return A list containing:
  #' \describe{
  #'  \item{ate}{The difference of weighted averages.}
  #'  \item{conf.int.pm}{The dispersion of the confidence interval (CI), wher CI = (ate - conf.int.pm, ate + conf.int.pm). }
  #' }
  #' 
  #' @examples
  #' data.example <- data.frame(V1 = rnorm(150), V2 = rnorm(150, 0.5, 1),
  #'  Treatment = c(rep(0, 100), rep(1, 50)),
  #'  Hospitalisation = c(rep(0,70), rep(1, 30), rep(1, 20), rep(0, 30)),
  #'  Death = (runif(150) < 0.1) * 1)
  #' 
  #' # # Match minority class with majority class with no trimming in the minority and 50% in the majority class
  #' matching.1 <- trimmedMatching(data.example, as.matrix(dist(data.example[, 1:2])),
  #'  trimm.level.majority = 0.5, strata = "Treatment", strata.levels = c(0, 1))
  #'  
  #' ate <- weightedMeansDifference(data.example$Hospitalisation[matching.1$trimmedMajorityClass],
  #'  data.example$Hospitalisation[matching.1$trimmedMinorityClass],
  #'   matching.1$trimmedMajorityClass.weight, matching.1$trimmedMinorityClass.weight)
  #'  
  #' @references H Inouzhe, I Barrio, MX Rodríguez-Álvarez, P Gordaliza, I Bengoechea and JM Quintana. (2023) A study on group fairness in healthcare outcomes for nursing home residents during the COVID-19 pandemic in the Basque Country.
  #'
  #' @export
  #'
  
  V1.X <- sum(weights.X)
  V2.X <- sum(weights.X^2)
  weighted_mean.X <- sum(weights.X * X) / V1.X
  V1.Y <- sum(weights.Y)
  V2.Y <- sum(weights.Y^2)
  weighted_mean.Y <- sum(weights.Y * Y) / V1.Y
  var.X <- var(X)/(V1.X^2 / V2.X)
  var.Y <- var(Y)/(V1.Y^2 / V2.Y)
  conf.int <- -qnorm(confidence.level/2) * sqrt(var.X + var.Y)
  return(list(ate = weighted_mean.Y - weighted_mean.X, conf.int.pm = conf.int))
}
# docstring::docstring(weightedMeansDifference)