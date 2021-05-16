#### EM Discret ####

rm(list = ls())


# Palette
library(viridisLite)
c.pal <- viridis(2)

# Forward/Backward Algorithm #

# Fonction qui renvoie la coordonnée de l'obersvation V
which_vec <- function(seq, V) {
  s <- vector(length = length(seq))
  for (i in 1:length(seq)) {
    s[i] <- which(seq[i] == V)
  }
  return(s)
}

# ----------------- Backward -----------------

backward_prob <- function(A, B, E, V, mu, seq) {
  s <- which_vec(seq, V)
  # Initialisation
  Beta_T <- rep(1, length(E))

  ans <- cbind(matrix(NA, nrow = length(E), ncol = length(s) - 1), Beta_T)

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

# On veut calculer alpha_t

# t représente la fin
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


# ---- Implementation des quantites necessaires aux formules de reestimation de Baum_Welch ----

# Fonction xi

xi <- function(A, B, E, V, mu, O, i, j, k) {
  return((forward_prob(A, B, E, V, mu, O)[i, k] *
            A[i, j] * B[j, which_vec(O[k + 1], V)] * 
            backward_prob(A, B, E, V, mu, O)[j, k + 1])
  / prob(0, A, B, E, V, mu, O))
}


# fonction gamma

gamma <- function(A, B, E, V, mu, O, i, k) {
  return((forward_prob(A, B, E, V, mu, O)[i, k] * backward_prob(A, B, E, V, mu, O)[i, k]) / prob(0, A, B, E, V, mu, O))
}

# fonction mu(i)
mu_estim_i <- function(A, B, E, V, mu, O, i, m) {
  mu_i <- 0
  for (l in 1:m) {
    mu_i <- mu_i + gamma(A, B, E, V, mu, O, i, 1) # pour plusieurs sequences O[,l]
  }
  return(mu_i / m)
}

# aij
a_ij_estim <- function(A, B, E, V, mu, O, i, j) {
  numerateur <- numeric(length = length(O) - 1)
  for (k in (1:length(O) - 1)) { # O[, l]
    numerateur[k] <- xi(A, B, E, V, mu, O, i, j, k) # O [, l]
  }

  numerateur <- sum(numerateur)

  denominateur <- numeric(length = length(O) - 1)
  for (k in 1:length(O) - 1) { # O[, l]
    denominateur[k] <- gamma(A, B, E, V, mu, O, i, k) # O[, l]
  }
  denominateur <- sum(denominateur)
  return(numerateur / denominateur)
}


b_jl_estim <- function(A, B, E, V, mu, O, i, v) {
  s <- which_vec(O, V)
  numerateur <- numeric(length = length(O) - 1)
  for (k in 2:length(O) - 1) {
    if (s[k] == v) { # indicatrice (la cordonnee k de la sequence O correspond aux symboles v in V)
      numerateur[k] <- gamma(A, B, E, V, mu, O, i, k)
    }
  }
  numerateur <- sum(numerateur)


  denominateur <- numeric(length = length(O) - 1)
  for (k in (2:length(O) - 1)) {
    denominateur[k] <- gamma(A, B, E, V, mu, O, i, k)
  }
  denominateur <- sum(denominateur)


  return(numerateur / denominateur)
}

# fonction qui retourne la matrice A avec les coef reestimes

A_estim <- function(A, B, E, V, mu, O) {
  A_new <- matrix(nrow = length(E), ncol = length(E))
  for (u in 1:length(E)) {
    for (v in 1:length(E)) {
      A_new [u, v] <- a_ij_estim(A, B, E, V, mu, O, u, v)
    }
  }
  return(A_new)
}

# fonction qui retourne la matrice B avec les coef reestimes

B_estim <- function(A, B, E, V, mu, O) {
  B_new <- matrix(nrow = length(E), ncol = length(V))
  for (u in 1:length(E)) {
    for (v in 1:length(V)) {
      B_new [u, v] <- b_jl_estim(A, B, E, V, mu, O, u, v)
    }
  }
  return(B_new)
}

# fonction qui renvoie le vecteur de proba initial reestime

mu_estim <- function(A, B, E, V, mu, O, m) {
  n <- length(E)
  mu_new <- numeric(length = n)
  for (u in 1:n) {
    mu_new[u] <- mu_estim_i(A, B, E, V, mu, O, u, m)
  }
  return(mu_new)
}

# ---- Algorithme de Baum Welch (cas discret) ----

Baum_Welch <- function(A0, B0, mu0, E, V, observation, max_iteration, eps) {
  prob0 <- prob(length(observation), A0, B0, E, V, mu0, observation)
  p <- 1
  mu1 <- mu_estim(A0, B0, E, V, mu0, observation, 1)
  A1 <- A_estim(A0, B0, E, V, mu0, observation)
  B1 <- B_estim(A0, B0, E, V, mu0, observation)
  prob1 <- prob(length(observation), A1, B1, E, V, mu1, observation)
  if (is.na(prob1) == TRUE) {
    return(list(proba_observation = prob1, proba_initial = mu1, A = A1, B = B1, nb_iteration = p))
  }
  while (abs(prob0 - prob1) > eps & p < max_iteration) {
    prob0 <- prob1
    mu0 <- mu1
    A0 <- A1
    B0 <- B1
    mu1 <- mu_estim(A0, B0, E, V, mu0, observation, 1)
    A1 <- A_estim(A0, B0, E, V, mu0, observation)
    B1 <- B_estim(A0, B0, E, V, mu0, observation)
    prob1 <- prob(length(observation), A1, B1, E, V, mu1, observation)
    if (is.na(prob1) == TRUE) {
      return(list(proba_observation = prob1, proba_initial = mu1, A = A1, B = B1, nb_iteration = p))
    }
    p <- p + 1
  }

  return(list(proba_observation = prob1, proba_initial = mu1, A = A1, B = B1, nb_iteration = p))
}


# ---- Application de l'algorithme de Baum-Welch ----

A <- matrix(c(
  0.6, 0.3, 0.1,
  0.1, 0.8, 0.1,
  0.1, 0.3, 0.6
),
byrow = TRUE, ncol = 3
)

B <- matrix(c(
  1, .5,
  0, 0,
  0.5, 1
), ncol = 2)

E <- c(1, 2, 3)
V <- c("a", "b")
mu <- c(0.6, 0.2, 0.2)
observation <- c("a", "b", "b", "a")

(proba_initiale <- prob(length(observation), A, B, E, V, mu, observation))

# Estimation des parametres

HMM_result <- Baum_Welch(A, B, mu, E, V, observation, 15, 1e-8)

# Evolution de la probabilite en fonction du nombre d iterations

list_proba <- function(A, B, mu, E, V, obsrvation, eps, ite) {
  proba <- numeric(length = ite)
  for (i in 1:ite) {
    proba[i] <- Baum_Welch(A, B, mu, E, V, observation, i, eps)$proba_observation
  }
  return(proba)
}

evol_proba <- list_proba(A, B, mu, E, V, observation, 1e-8, 15)
plot(evol_proba,
  type = "l", col = c.pal[1], ylab = "probabilité", xlab = "nombre d'itérations",
  main = "Evolution de la probabilité d'observation"
)
