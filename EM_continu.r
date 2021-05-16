# EM cas gaussien  #

rm(list = ls())

# install.package("HiddenMarkov")
library(HiddenMarkov)

# Palette
library(viridisLite)
c.pal <- viridis(2)

set.seed(666)

# Probabilité forward-backward cas gaussien

proba <- function(A, mu, moyenne, sigma, sequence) {
  y <- forwardback(sequence, Pi = A, mu, "norm", list(mean = moyenne, sd = sigma))
  alpha <- exp(y$logalpha)
  beta <- exp(y$logbeta)
  return(list(alpha = alpha, beta = beta))
}

# Initialisation des parametres du HMM

A <- matrix(c(
  0.6, 0.3, 0.1, 
  0.1, 0.8, 0.1,
  0.1, 0.3, 0.6
),
byrow = TRUE, nrow = 3
)

mu <- c(0.2, 0.6, 0.2)

x <- dthmm(
  NULL, A, mu, "norm",
  list(mean = c((-2), 0, 2), sd = c(0.5, 0.5, 0.5))
)


# Simulation d'une sequence de longueur 1000 issu des parametres consideres

x <- simulate(x, nsim = 1000)

# Application de Baum_welch sur la sequence o simulee precedemment

bw_control <- bwcontrol(
  maxiter = 30, tol = 1e-08, prt = TRUE, posdiff = TRUE,
  converge = expression(diff < tol)
)

y <- BaumWelch(x, control = bw_control)

# Evolution de l'estimation des parametres en fonction des iterations

parameter_list <- function(A, mu, n) {
  m_1 <- numeric(length = n)
  m_2 <- numeric(length = n)
  m_3 <- numeric(length = n)

  x <- dthmm(
    NULL, A, mu, "norm",
    list(mean = c((-2), 0, 2), sd = c(0.5, 0.5, 0.5))
  )
  x <- simulate(x, nsim = 1000)

  for (i in 1:n) {
    bw_control <- bwcontrol(
      maxiter = i, tol = 1e-15, prt = FALSE,
      posdiff = FALSE, converge = expression(diff < tol)
    )

    estim <- BaumWelch(x, bw_control)
    m_1[i] <- estim$pm$mean[1]
    m_2[i] <- estim$pm$mean[2]
    m_3[i] <- estim$pm$mean[3]
  }
  return(list(m_1 = m_1, m_2 = m_2, m_3 = m_3))
}


parameter <- parameter_list(A, mu, 30)

par(mfrow = c(3, 1))

plot(parameter$m_1,
  type = "l", col = c.pal[1],
  xlab = "nombre d'itérations",
  main = expression(paste("estimation de ", m[1]))
)

plot(parameter$m_2,
  type = "l", col = c.pal[1],
  xlab = "nombre d'itérations",
  main = expression(paste("estimation de ", m[2]))
)

plot(parameter$m_3,
  type = "l", col = c.pal[1],
  xlab = "nombre d'itérations",
  main = expression(paste("estimation de ", m[3]))
)


# Evolution de la log-vraisemblance complete en fonction des itetations

par(mfrow = c(1, 1))

likelihood_list <- function(A, mu, n) {
  likelihood <- numeric(length = n)
  x <- dthmm(
    NULL, A, mu, "norm",
    list(mean = c((-2), 0, 2), sd = c(0.5, 0.5, 0.5))
  )
  x <- simulate(x, nsim = 1000)

  for (i in 1:n) {
    bw_control <- bwcontrol(maxiter = i, tol = 1e-15, prt = FALSE, posdiff = FALSE, converge = expression(diff < tol))

    estim <- BaumWelch(x, bw_control)
    likelihood[i] <- logLik(estim)
  }
  return(likelihood)
}


likelihood <- likelihood_list(A, mu, 30)

plot(likelihood, type = "l", col = c.pal[1], xlab = "nombre d'itérations", main = "Evolution de la log-vraisemblance complète")


# Construction Bootstrap IC

# Première estimation des paramètres afin de générer K séquences :

A <- matrix(c(
  0.6, 0.3, 0.1, 0.1, 0.8, 0.1,
  0.1, 0.3, 0.6
),
byrow = TRUE, nrow = 3
)

mu <- c(0.2, 0.6, 0.2)

x <- dthmm(
  NULL, A, mu, "norm",
  list(mean = c((-2), 0, 2), sd = c(0.5, 0.5, 0.5))
)

x <- simulate(x, nsim = 1000)
bw_control <- bwcontrol(maxiter = 30, tol = 1e-08, prt = FALSE, posdiff = TRUE, converge = expression(diff < tol))
y <- BaumWelch(x, control = bw_control)

bootstrap_1 <- dthmm(NULL, Pi = y$Pi, delta = y$delta, "norm", list(mean = c((-1.96), (-0.04), 1.96), sd = c(0.48, 0.49, 0.5)))

# On simule K sequences et on effectue K fois EM

bootstrap_IC <- function(HMM, K) {
  parameter_list2 <- matrix(NA, ncol = 2, nrow = K)
  MSE <- matrix(NA, ncol = 2 , nrow = K)

  for (i in 1:K) {
    HMM <- simulate(HMM, nsim = 1000)

    bw_control <- bwcontrol(maxiter = 30, tol = 1e-08, prt = FALSE, posdiff = TRUE, converge = expression(diff < tol))

    HMM_new <- BaumWelch(HMM, control = bw_control)

    parameter_list2[i, 1] <- HMM_new$pm$mean[3]
    parameter_list2[i, 2] <- HMM_new$Pi[2, 2]
    
    MSE[i,1] <- mean((parameter_list2[1:i,1]-rep(2,i))^2)
    MSE[i,2] <- mean((parameter_list2[1:i,2]-rep(.8,i))^2)
  }
  return(list(parameter_list2 = parameter_list2,MSE = MSE))
}

echantillon_bootstrap <- bootstrap_IC(bootstrap_1, 500)

hist(echantillon_bootstrap$parameter_list2[, 1],
  breaks = 30, col = c.pal[1], border = c.pal[2], ylab = "Fréquence", xlab = expression(paste("Estimation de ", m[3])),
  main = expression(paste("Histogramme de l'echantillon bootstrap de ", m[3]))
)

hist(echantillon_bootstrap$parameter_list2[, 2],
  breaks = 30, col = c.pal[1], border = c.pal[2], ylab = "Fréquence", xlab = expression(paste("Estimation de ", a[22])),
  main = expression(paste("Histogramme de l'echantillon bootstrap de ", a[22]))
)

# Calcul de l'intervalle de confiance avec les quantiles empiriques

born_inf <- sort(echantillon_bootstrap$parameter_list, na.last = NA)[floor((0.05 * 500) / 2) + 1]
born_sup <- sort(echantillon_bootstrap$parameter_list, na.last = NA)[floor((1 - 0.05 / 2) * 500)]

# MSE

plot(echantillon_bootstrap$MSE[,1],
     col = c.pal[1], border = c.pal[2], ylab = "MSE", xlab = "Itéaration",
     main = expression(paste("Histogramme de l'echantillon bootstrap de ", m[3]))
)

plot(echantillon_bootstrap$MSE[,1], 
     type = "l",ylab = "MSE", 
     xlab = "Itération",
     col = c.pal[1],
     main =expression(paste("Histogramme de l'echantillon bootstrap de ", m[3])))


plot(echantillon_bootstrap$MSE[,2], 
     type = "l",ylab = "MSE", 
     xlab = "Itération",
     col = c.pal[1],
     main =expression(paste("Histogramme de l'echantillon bootstrap de ", a[22])))


