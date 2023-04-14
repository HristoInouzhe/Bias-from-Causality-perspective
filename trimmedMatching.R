trimmedMatching <- function(data, distance_matrix, subset_to_use = c(),
                            trimm.level.majority = 0,  trimm.level.minority = 0,
                            strata, strata.levels){
  #' trimmedMatching
  #'
  #' Computes a trimmed optimal transport match between two sets of points and returns a pair of weighted subsamples based on the matching.
  #'
  #' @param data A data frame, with data entries as rows, containing both sets of data points.
  #' @param distance_matrix The cost matrix, usually a power of the distance matrix, between the entries of data.
  #' @param subset_to_use A subset of data for which we want to make the matching. Default is the whole data set.
  #' @param trimm.level.majority Trimming level to use for the majority class (set), must be strictly greater than 0 and less than 1.
  #' @param trimm.level.minority Trimming level to use for the minority class (set), must be strictly greater than 0 and less than 1.
  #' @param strata The name of the dicotomic variable, as in the data frame data, which defines the two classes to be matched.
  #' @param strata.levels The two levels of the dicotomic variable strata.
  #' 
  #' @return A list containing:
  #' \describe{
  #'  \item{trimmedMajorityClass}{Indexes of the points in the data frame data of the majority class that are left after the trimmed optimal matching.}
  #'  \item{trimmedMajorityClass.weight}{The weights associated with the majority class left after the matching. Positive and sum to one.}
  #'  \item{trimmedMinorityClass}{Indexes of the points in the data frame data of the minority class that are left after the trimmed optimal matching.}
  #'  \item{trimmedMinorityClass.weight}{The weights associated with the minority class left after the matching. Positive and sum to one.}
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
  #' # # Match minority class with majority class with no trimming in the minority and 33.33% in the majority class for hospitalisations.
  #' matching.2 <- trimmedMatching(data.example, as.matrix(dist(data.example[, 1:2])), subset_to_use = which(data.example$Hospitalisation == 1),
  #'  trimm.level.majority = 1/3, strata = "Treatment", strata.levels = c(0, 1))
  #'  
  #' @references H Inouzhe, I Barrio, MX Rodríguez-Álvarez, P Gordaliza, I Bengoechea and JM Quintana. (2023) A study on group fairness in healthcare outcomes for nursing home residents during the COVID-19 pandemic in the Basque Country.
  #'
  #' @export
  #'
  
  if(length(subset_to_use) == 0){
    subset_to_use <- 1:(dim(data)[1])
  }
  data.majority.class <- subset_to_use[data[[strata]][subset_to_use] == strata.levels[1]]
  data.minority.class <- subset_to_use[data[[strata]][subset_to_use] == strata.levels[2]]
  N_0 <- length(data.majority.class)
  N_1 <- length(data.minority.class)
  
  alpha_1 <- trimm.level.majority
  alpha_2 <- trimm.level.minority
  
  P <- c(rep(1 / ((1 - alpha_1) * N_0), N_0), alpha_2 / (1 - alpha_2))
  Q <- c(rep(1 / ((1 - alpha_2) * N_1), N_1), alpha_1 / (1 - alpha_1))
  
  cost.matrix <- distance_matrix[data.majority.class, data.minority.class]
  
  cost.matrix <- cbind(cost.matrix, rep(0, N_0))
  cost.matrix <- rbind(cost.matrix,
                       c(rep(0, N_1), ((1 - alpha_1) / (alpha_1 * (1 - alpha_2))) * sum(cost.matrix)))
  
  # t0 <- Sys.time()
  trimmed.transport <- transport::transport(P, Q, cost.matrix)
  # t1 <- Sys.time()
  # print(t1 - t0)
  
  trimmed.transport_1 <- trimmed.transport[trimmed.transport$from != (N_0 + 1), ]
  trimmed.transport_1 <- trimmed.transport_1[trimmed.transport_1$to != (N_1 + 1), ]
  
  non.trimmed.0 <- unique(trimmed.transport_1$from)
  non.trimmed.1 <- unique(trimmed.transport_1$to)
  
  N_0_1 <- length(non.trimmed.0)
  N_1_1 <- length(non.trimmed.1)
  
  origin <- data.majority.class[non.trimmed.0]
  destination <- data.minority.class[non.trimmed.1]
  
  trimmed.transport_1$from <- factor(trimmed.transport_1$from, levels = non.trimmed.0)
  trimmed.transport_1$to <- factor(trimmed.transport_1$to, levels = non.trimmed.1)
  
  origin.weight <- (trimmed.transport_1 %>% group_by(from) %>% summarise(weight = sum(mass)))$weight
  origin.weight <- origin.weight / sum(origin.weight)
  
  destination.weight <- (trimmed.transport_1 %>% group_by(to) %>% summarise(weight = sum(mass)))$weight
  destination.weight <- destination.weight / sum(destination.weight)
  
  return(list(trimmedMajorityClass = origin,
              trimmedMajorityClass.weight = origin.weight,
              trimmedMinorityClass = destination,
              trimmedMinorityClass.weight = destination.weight))
}
# docstring::docstring(trimmedMatching)