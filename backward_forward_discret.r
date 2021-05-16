rm(list = ls())

# Initialisation
A <- matrix(c(
  0.3, 0.5, 0.2,
  0, 0.3, 0.7,
  0, 0, 1
), byrow = TRUE, ncol = 3)

B <- matrix(c(
  1, 0,
  0.5, 0.5,
  0, 1
), byrow = TRUE, ncol = 2)

mu <- c(0.6, 0.4, 0)
E <- c(1, 2, 3)
V <- c("a", "b")
seq <- c("a", "b", "a", "a", "b")

###### --------------------------------- CAS DISCRET --------------------------------------------

# Fonction qui renvoie la coordonnée de l'obersvation V
which_vec <- function(seq, V) {
  s <- vector(length = length(seq))
  for (i in 1:length(seq)) {
    s[i] <- which(seq[i] == V)
  }
  return(s)
}


# -------------- Backward --------------
backward_prob <- function(A, B, E, V, mu, seq) {
  s <- which_vec(seq, V)

  # Initialisation
  Beta_1 <- rep(1, length(E))
  ans <- cbind(matrix(NA, nrow = length(E), ncol = length(s) - 1), Beta_1)

  # Boucle sur le temps ``
  l <- length(s) - 1
  while (l >= 1) {
    ans[, l] <- A %*% (B[, s[l + 1]] * ans[, l + 1])
    l <- l - 1
  }
  # Resultat au temps t=0 final (utilisé plus tard)
  res <- (mu * B[, s[1]]) %*% ans[, 1]

  return(ans)
}


# ----------------- Forward -----------------
forward_prob <- function(A, B, E, V, mu, seq) {
  s <- which_vec(seq, V)

  # Initialisation
  alpha_1 <- mu * B[, s[1]]

  ans <- cbind(alpha_1, matrix(NA, nrow = length(E), ncol = length(s) - 1))

  # Boucle sur le temps
  for (l in 2:length(s) - 1) {
    ans[, l + 1] <- (t(A) %*% ans[, l]) * B[, s[l + 1]]
  }
  return(ans)
}


# ------------ Proba au temps t en utilisant backward-forward ---------
prob <- function(t, A, B, E, V, mu, seq) {
  if (t == 0) {
    s <- which_vec(seq, V)
    beta_t <- backward_prob(A, B, E, V, mu, seq)
    return((mu * B[, s[1]]) %*% beta_t[, 1])
  }
  else {
    alpha_t <- forward_prob(A, B, E, V, mu, seq)
    beta_t <- backward_prob(A, B, E, V, mu, seq)
    return(sum(alpha_t[, t] * beta_t[, t]))
  }
}


(alpha_t <- forward_prob(A, B, E, V, mu, seq))
(beta_t <- backward_prob(A, B, E, V, mu, seq))

prob(2, A, B, E, V, mu, seq)
