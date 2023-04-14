prob_outcome_weighted <- function(X, Y, weights.X, weights.Y, test.level = 0.05){

  V1.X <- sum(weights.X)
  V2.X <- sum(weights.X^2)
  weighted_mean.X <- (weights.X %*% X)[1,1] / V1.X
  V1.Y <- sum(weights.Y)
  V2.Y <- sum(weights.Y^2)
  weighted_mean.Y <- (weights.Y %*% Y)[1,1] / V1.Y
  var.X <- var(X)/(V1.X^2 / V2.X)
  var.Y <- var(Y)/(V1.Y^2 / V2.Y)
  conf.int <- -qnorm(test.level/2) * sqrt(var.X + var.Y)
  return(list(ate = weighted_mean.Y - weighted_mean.X, conf.int.pm = conf.int))
}