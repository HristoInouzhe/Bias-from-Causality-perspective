wasserDistance <- function(p1, p2, cost.matrix, d){
  #' wasserDistance
  #'
  #' Computes the d-Wasserstein distance.
  #'
  #' @param p1 A 1D array of probabilities.
  #' @param p2 A 1D array of probabilities.
  #' @param cost_matrix The d-power of the Euclidean distance between two sets of points corresponding with p1 and p2, respectively..
  #' @param d The power d for the d-Wasserstein distace. 
  #' 
  #' @return An estimate of the d-Wasserstein distance.  
  #' 
  #' @examples
  #' data.example <- data.frame(V1 = rnorm(150), V2 = rnorm(150, 0.5, 1),
  #'  Treatment = c(rep(0, 100), rep(1, 50)),
  #'  Hospitalisation = c(rep(0,70), rep(1, 30), rep(1, 20), rep(0, 30)),
  #'  Death = (runif(150) < 0.1) * 1)
  #' 
  #' 2.Wasser.Dist <- wasserDistance(rep(1/100, 100), rep(1/50, 50), (as.matrix(dist(data.example[, 1:2]))^2)[1:100, 101:150], 2)
  #'  
  #' @references H Inouzhe, I Barrio, MX Rodríguez-Álvarez, P Gordaliza, I Bengoechea and JM Quintana. (2023) A study on group fairness in healthcare outcomes for nursing home residents during the COVID-19 pandemic in the Basque Country.
  #'
  #' @export
  #'
  
  trsp.ot <- transport::transport(p1, p2, cost.matrix)
  wdist <-((cost.matrix[cbind(trsp.ot$from, trsp.ot$to)]%*%trsp.ot$mass)[1,1])^(1/d)
  return(wdist)
}
# docstring::docstring(wasserDistance)