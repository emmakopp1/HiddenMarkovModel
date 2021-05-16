rm(list = ls())

install.packages("HiddenMarkov")
library(HiddenMarkov)

# Paramètres
E <- c(1, 2, 3)

A <- matrix(c(
  0.6, 0.3, 0.1, 
  0.1, 0.8, 0.1,
  0.1, 0.3, 0.6
),
byrow = TRUE, nrow = 3
)

mu <- c(.2, 0.6, .2)
B <- c(-2, 0, 2)
sigma <- c(1, 1, 1)
n <- 5


# Backward - Forward
backwardforward <- function(A, mu, moyenne, sigma, sequence) {
  y <- forwardback(sequence, Pi = A, mu, "norm", list(mean = moyenne, sd = sigma))

  alpha <- t(exp(y$logalpha)) # le package renvoie la transposée de la matrice souhaitée
  beta <- t(exp(y$logbeta))
  return(list(alpha = alpha, beta = beta))
}

# Application
sequence <- c(1.73, 2.17, 1.28, 1.44, 1.51)

sequence <- sample(seq(-10, 10, length.out = 10), size = 1000,replace = TRUE)
alpha <- backwardforward(A, mu, B, sigma, sequence)$alpha
beta <- backwardforward(A, mu, B, sigma, sequence)$beta


round(alpha, digits = 4)
round(beta, digits = 4)

# Probabilité d'observation de la séquence :
t <- 3
print(alpha[, t] %*% beta[, t])
