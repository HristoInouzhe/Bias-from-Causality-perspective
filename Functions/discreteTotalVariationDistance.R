discreteTotalVariationDistance <- function(x, y, w_x = NULL, w_y = NULL){
  #' discreteTotalVariationDistance
  #'
  #' Computes a 1D Total Variation distances (TVd) estimate between two classes for discrete (categorival) data.
  #'
  #' @param x A 1D array.
  #' @param y A 1D array.
  #' @param w_x Weights for x, will be normalised to one.
  #' @param w_y Weights for x, will be normalised to one. 
  #' 
  #' @return An estimate of the total variation distance.  
  #' 
  #' @examples
  #' data.example <- data.frame(V1 = rnorm(150), V2 = rnorm(150, 0.5, 1),
  #'  Treatment = c(rep(0, 100), rep(1, 50)),
  #'  Hospitalisation = c(rep(0,70), rep(1, 30), rep(1, 20), rep(0, 30)),
  #'  Death = (runif(150) < 0.1) * 1)
  #' 
  #' TVd <- discreteTotalVariationDistance(data.example$Death[data.example$Treatment == 0],
  #'  data.example$Death[data.example$Treatment == 1])
  #'  
  #' @references H Inouzhe, I Barrio, MX Rodríguez-Álvarez, P Gordaliza, I Bengoechea and JM Quintana. (2023) A study on group fairness in healthcare outcomes for nursing home residents during the COVID-19 pandemic in the Basque Country.
  #'
  #' @export
  #'

  if(is.null(w_x)){
    w_x <- rep(1, length(x))
  } 
  if(is.null(w_y)){
    w_y <- rep(1, length(y))
  } 
  probability.1.0 <- data.frame(V0 = x, V1 = w_x)
  probability.2.0 <- data.frame(V0 = y, V1 = w_y)
  probability.1 <- probability.1.0 %>% group_by(V0) %>% summarise(sum = sum(V1) / sum(probability.1.0$V1))
  probability.2 <- probability.2.0 %>% group_by(V0) %>% summarise(sum = sum(V1) / sum(probability.2.0$V1))
  
  probability.join <- left_join(probability.1, probability.2, by = "V0")
  probability.join$sum.x[is.na(probability.join$sum.x)] <- rep(0, sum(is.na(probability.join$sum.x)))
  probability.join$sum.y[is.na(probability.join$sum.y)] <- rep(0, sum(is.na(probability.join$sum.y)))
  
  return(
    0.5 * sum(abs(probability.join$sum.x - probability.join$sum.y))
  )
}
docstring::docstring(discreteTotalVariationDistance)