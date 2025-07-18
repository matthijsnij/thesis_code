compute_z1 <- function(x1, x2, eps) {
  ifelse(x1 <= 0.5,
         5 * sin(2 * pi * x1) + 3 * x2 + eps,
         5 * sin(2 * pi * x1) - 3 * x2 + eps)
}

compute_z2 <- function(x1, x2, eps) {
  ifelse(x2 <= 0.3,
         4 * cos(pi * x2) + x1 + eps,
         8 * cos(pi * x2) - x1 - 2 + eps)
}

set.seed(42)  # for reproducibility

# Number of observations
n <- 10000

# Generate uniform features and noise
x1 <- runif(n)
x2 <- runif(n)
eps1 <- rnorm(n)
eps2 <- rnorm(n)

z1 <- compute_z1(x1, x2, eps1)
z2 <- compute_z2(x1, x2, eps2)
z <- cbind(z1,z2)

assign_class <- function(z_matrix) {
  apply(z_matrix, 1, function(z) {
    if (all(z < 0)) {
      return(0)
    } else if (z[1] > z[2]) {
      return(1)
    } else {
      return(2)
    }
  })
}

y <- assign_class(z)
table(y) / length(y)  # Proportion of each class