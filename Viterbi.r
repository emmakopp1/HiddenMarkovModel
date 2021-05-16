##### Viterbi ####

rm(list = ls())

# Initialisation des parametres du HMM

A <- matrix(c(0.3, 0, 0, 1 / 2, 0.3, 0, 0.2, 0.7, 1), ncol = 3) # matrice de transition des etats caches
B <- matrix(c(1, 1 / 2, 0, 0, 0.5, 1), ncol = 2) # matrice d'emission des etats observes
mu <- c(0.6, 0.4, 0) # vecteur de proba initiale
E <- c(1, 2, 3) # ensemble des etats caches qu'on code 1,2,3
V <- c("a", "b") # ensemble des etats observables


# sequence des observations

o <- c("a", "a", "b", "b")

# fonction permettant de transformer des observations catégoriques en numeriques
which_vec <- function(seq, V) {
  s <- vector(length = length(seq))
  for (i in 1:length(seq)) {
    s[i] <- which(seq[i] == V)
  }
  return(s)
}

# Algorithme de Viterbi applique à la sequence o

viterbi <- function(A, B, E, V, mu, s) {
  s <- which_vec(s, V)
  delta <- mu * B [, s[1]]
  phi <- rep(0, length(E))
  delta_matrix <- cbind(delta, matrix(NA, nrow = length(E), ncol = length(s) - 1))
  phi_matrix <- cbind(phi, matrix(NA, nrow = length(E), ncol = length(s) - 1))
  t <- 2
  while (t <= length(s)) {
    j <- 1
    while (j <= length(E)) {
      delta_matrix[j, t] <- max(delta_matrix[, t - 1] * A[, j]) * B[j, s[t]]
      phi_matrix[j, t] <- which.max(delta_matrix[, t - 1] * A[, j])
      j <- j + 1
    }
    t <- t + 1
  }
  print(delta_matrix)
  print(phi_matrix)
  prob_max <- max(delta_matrix[, length(s)])
  hidden_list <- vector(length = length(s))
  hidden_list[length(s)] <- which.max(delta_matrix[, length(s)])
  t <- length(s) - 1
  while (t >= 1) {
    print(t)
    hidden_list[t] <- phi_matrix[hidden_list[t + 1], t + 1]
    t <- t - 1
  }
  return(list(proba_max = prob_max, hidden_list = hidden_list))
}

viterbi <- viterbi(A, B, E, V, mu, o)

chemin_optimal <- viterbi$hidden_list
proba_optimal <- viterbi$proba_max
