#### Gibbs ####

rm(list = ls())
set.seed(66)

#install.packages("DirichletReg")
library(DirichletReg)
library(HiddenMarkov)



## Simulation loi conditionnelle complète pour proba initiale mu

mu_sim_cond <- function(E, h) {
  vect <- rep(1, times = length(E))
  for (i in 1:length(E)) {
    vect[i] <- vect[i] + as.numeric(E[i] == h[1])
  }
  return(rdirichlet(1, vect))
}


# n_ij : nb de transition entre i et j
n <- function(E, i, j, h) {
  s <- 0
  for (u in 2:length(h)) {
    s <- s + as.numeric(((h[u - 1] == E[i]) & (h[u] == E[j])))
  }
  return(s)
}

## Simulation loi conditionnelle pour ligne de A
row_sim_cond <- function(E, h, i) {
  vect <- rep(1, times = length(E))
  for (u in 1:length(E)) {
    vect[u] <- vect[u] + n(E, i, u, h)
  }
  return(rdirichlet(1, vect))
}


# Estimation de A 
A_estim_cond <- function(E, h) {
  A <- matrix(data = NA, ncol = length(E), nrow = length(E))
  for (i in 1:length(E)) {
    A[i, ] <- row_sim_cond(E, h, i)
  }
  return(A)
}


# simulation de la loi conditionnelle complete des mi

# Si : moyenne des observations emise par l'état e_j
S <- function(E, o, h, i) {
  s <- 0
  for (u in 1:length(h)) {
    s <- s + as.numeric((h[u] == E[i])) * o[u]
  }
  return(s)
}


# n_i : nb de passage par l'etat ei
n_i <- function(E, i, h) {
  s <- 0
  for (u in 1:length(h)) {
    s <- s + as.numeric(((h[u] == E[i])))
  }
  return(s)
}



mi_vect <- function(E, h, o, sigma) {
  u <- (max(o) + min(o)) / 2
  w <- (max(o) - min(o))^(-2)
  
  vect <- numeric(length = length(E))
  
  for (i in 1:length(E)) {
    # moyenne estimée 
    m <- (S(E, o, h, i) + w * u * sigma) / (n_i(E, i, h) + w * sigma)
    # variance estimée 
    v <- sigma / (n_i(E, i, h) + w * sigma)
    vect[i] <- rnorm(1, m, sqrt(v))
  }
  return(vect)
}


# simulation de la loi conditionnelle complete de sigma^(-2)
sigma_cond <- function(E, h, o, sigma, alpha, beta,moyenne) {
  v <- 0
  for (k in 1:length(o)) {
    v <- v + (o[k] - moyenne[h[k]])^2
  }
  m <- alpha + length(h) / 2
  v <- (1 / 2) * v + beta
  return(rgamma(1, m, v))
}


# Simulation de la loi conditionnelle complete de beta
beta_sim_cond <- function(o, sigma_inv, f =.2, g =10 / (max(o) - min(o))^2, alpha =2) {
  return(rgamma(1, f + alpha, g + sigma_inv))
}

#beta_sim_cond(o, sigma_cond(E, h, o, 0.3, 2, beta = 4))

## Fonction forward / backward continue

# Probabilité forward backward cas continu
proba_2 <- function(A, mu, moyenne, sigma, sequence) {
  y <- forwardback(sequence, Pi = A, mu, "norm", list(mean = moyenne, sd = sigma))
  alpha <- (y$logalpha)
  beta <- (y$logbeta)
  return(list(alpha = alpha, beta = beta))
}

## Simulation de la chaine de Markov non homogene

#Simulation via log proba
gumbel_sample <- function(a){
  g <- -log(-log(runif(3)))
  return(which.max(a + g))
}


# proba initiale
gene_markov_chain <- function(E, A, m_vect, mu, sigma, o) {
  markov_chain <- numeric(length = length(o))
  proba_ini <- numeric(length = length(E))
  for (i in 1:length(E)) { 
    proba_ini[i] <- log(mu[i]) + log(dnorm(o[1], mean = m_vect[i], sd = sqrt(sigma)))
    + proba_2(A, mu, m_vect, rep(sigma, times = length(E)), o)$beta[2, 1]
    
  }
  
  markov_chain[1] <- gumbel_sample(proba_ini)
  for (k in 2:length(o)) { 
    A_k <- matrix(data = NA, nrow = length(E), ncol = length(E))
    for (i in 1:length(E)) {
      for (j in 1:length(E)) {
        A_k[i, j] <- log(A[i, j]) +  dnorm(o[k], m_vect[j], sd = sqrt(sigma),log = TRUE) +
          proba_2(A, mu, m_vect, rep(sigma, times = length(E)), o)$beta[k, j]
      }
    }
    markov_chain[k] <- gumbel_sample(A_k[markov_chain[k - 1],])
  }
  return(markov_chain)
}


### Echantillonnage de Gibbs
Gibbs_sampler <- function(E, o, h, alpha, beta, sigma, ite) {
  parameter_test <- matrix(NA, ncol = 1, nrow = ite)
  
  o_matrix <- matrix(data = NA, nrow = length(o), ncol = ite)
  m_vect <- mi_vect(E, h, o, sigma)
  sigma_inv <- sigma_cond(E, h, o, sigma, alpha, beta,m_vect)
  beta <- beta_sim_cond(o, sigma_inv)
  A <- A_estim_cond(E, h)
  mu <- mu_sim_cond(E, h)
  markov_chain <- gene_markov_chain(E, A, m_vect, mu, sigma_inv^(-1), o)
  o_matrix[, 1] <- markov_chain
  parameter_test[1,] <- m_vect[2]
  
  for (i in 2:ite) {
    m_vect <- mi_vect(E, markov_chain,o, sigma_inv^(-1))
    sigma_inv <- sigma_cond(E, markov_chain, o, sigma_inv^(-1), alpha, beta, m_vect)
    beta <- beta_sim_cond(markov_chain, sigma_inv)
    A <- A_estim_cond(E, h)
    mu <- mu_sim_cond(E, h)
    markov_chain <- gene_markov_chain(E, A, m_vect, mu, sigma_inv^(-1), markov_chain)
    o_matrix[, i] <- markov_chain
    parameter_test[i,] <- m_vect[2]
  }
  return(list(simu_result = o_matrix, transition_matrix = A, initial_prob = mu, m_vect = m_vect, parameter = parameter_test))
}


# Exemple article 
E_article <- c(1, 2, 3)
A_article <- matrix(c(
  .6, .3, .1,
  0.1, .8, .1,
  0.1, 0.3, 0.6
),
byrow = TRUE, nrow = 3
)

mu_article <- c(.2, .6, .2)
B_article <- c(-2,0,2)
sigma_article <- .5
n_simu_article <- 1000


# Simulation manuelle de o et h selon les parametres A,B,mu
simulation <- function(E,n,mu,sigma,B,A){
  simu_hidden <- numeric(length = n)
  simu_observed <- numeric(length = n)
  
  simu_hidden <- sample(x = E, size = 1, prob = mu)
  simu_observed <- rnorm(1, mean = B[simu_hidden[1]], sd = sigma)
  
  for (i in 2:n) {
    simu_hidden[i] <- sample(x = E, size = 1, prob = A[simu_hidden[i - 1], ])
    simu_observed[i] <- rnorm(1, mean = B[simu_hidden[i]], sd = sigma)
  }
  return(list(simu_hidden = simu_hidden, simu_observed = simu_observed))
}

# Simulation HMM
test <- simulation(E_article,n_simu_article,mu_article,sigma_article,B_article,A_article)
simu_obs <- test$simu_observed
simu_hid <- test$simu_hidden

# Paramètres (c.f article)
# Moyenne initiale
m_vect_ini <- function(E,o){
  m_vect <- numeric(length = length(E))
  R <- max(o) - min(o)
  for (i in 1:length(E)){
    m_vect[i] <- min(o) + R/(length(E)) + (i - 1)*R/(length(E))
  }
  return(m_vect)
}

# Séquence cachée initiale
hidden_ini <- function(E,o,m_vect){
  hidden <- numeric(length = length(o))
  for (i in 1:length(o)){
    hidden[i] <- which.min((o[i]-m_vect)^2)
  }
  return(hidden)
}

# Variance initiale
sigma_ini <- function(o,h,m_vect){
  for (i in 1:length(o)){
    h[i] <- m_vect[h[i]]
  }
  return(mean((o - h)^2))
}

# Fonction d'initialisation des paramètres 
parameter_ini <- function(E,o){
  m_vect <- m_vect_ini(E,o)
  hidden <- hidden_ini(E,o,m_vect)
  sigma <- sigma_ini(o,hidden,m_vect)
  return(list(m_vect = m_vect,h = hidden,sigma = sigma))
}

# Initialisation des paramètres 
sigma_estim <- parameter_ini(E_article,simu_obs)$sigma
hidden_estim <- parameter_ini(E_article,simu_obs)$h
f <- 0.2
g <- 10 / (max(simu_obs) - min(simu_obs))^2
alpha <- 2
beta <- f/g
sigma <- 1 / rgamma(1, alpha, beta)

# Inference des parameters via Algorithm de Gibbs
par(mfrow = c(1,1))
gibbs_result <- Gibbs_sampler(E_article, simu_obs, hidden_estim, alpha, beta, sigma_article, ite = 2000)
plot(gibbs_result$parameter[,1], type = "l",ylim = c(-.5,.5), ylab = "Valeur" , xlab = "Itération",
     col = c.pal[1], main = expression(paste("Simulation de ",m[2])))

# Palette
library(viridisLite)
c.pal <- viridis(2)

hist(gibbs_result$parameter[500:2000,2], 
     breaks = 30, xlim = c(-.5,.5), 
     main = expression(paste("Histogramme des simulations de ",m[2])), 
     ylab = "Fréquence" , xlab = "Estimation du paramètree",
     col = c.pal[1], border = c.pal[2])
