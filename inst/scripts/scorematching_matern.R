library(rSPDE)

#alpha=1/range)1/kappa
#nu
set.seed(123)
nobs <- 101
x <- seq(from = 0, to = 1, length.out = nobs)
fem <- rSPDE.fem1d(x)
kappa <- 40
sigma <- 1
d <- 1
nu <- 2.6
tau <- sqrt(gamma(nu) / (kappa^(2 * nu) * (4 * pi)^(d / 2) *
                           gamma(nu + d / 2)))
range <- sqrt(8*nu)/kappa
op_cov <- matern.operators(
  loc_mesh = x, nu = nu, range = range, sigma = sigma,
  d = 1, m = 2, compute_higher_order = TRUE,
  parameterization = "matern"
)
v <- t(rSPDE.A1d(x, 0.5)) #Get the observation matrix for the point 0.5
c.true <- matern.covariance(abs(x - 0.5), kappa, nu, sigma) #covariance between each point and 0.5
Q <- rspde.matern.precision(
  kappa = kappa, nu = nu, tau = tau, rspde.order = 2, d = 1,
  fem_mesh_matrices = op_cov$fem_mesh_matrices
) #precision matrix, why is dimension of Q 3nx3n? -> (order+1)*nobs
A <- Diagonal(nobs)
Abar <- cbind(A, A, A)
w <- rbind(v, v, v)
c.approx_cov <- (Abar) %*% solve(Q, w)

# plot the result and compare with the true Matern covariance
plot(x, matern.covariance(abs(x - 0.5), kappa, nu, sigma),
     type = "l", ylab = "C(h)",
     xlab = "h", main = "Matern covariance and rational approximations"
)
lines(x, c.approx_cov, col = 2)



model <- rspde.matern(mesh=x,nu=nu,)


############################# INLA ######################################


# INLA code
inla.seed = sample.int(n=1E6, size=1)

sim_loc = matrix(c(0,0,3,3, 0, 3, 3, 0), nrow = 4, byrow = T)
mesh_sim = inla.mesh.2d(loc = sim_loc, max.edge=c(0.5, 1))
plot(mesh_sim)
mesh_sim$n

spde = inla.spde2.matern(mesh_sim, alpha = 2)

Q = inla.spde.precision(spde, theta=spde$param.inla$theta.initial)
omega_s = inla.qsample(n=1, Q = Q, seed = inla.seed)
omega_s = omega_s[ ,1]  # Spatial random field

# Simulated data
n = 1000
loc = matrix(runif(2*n), n)*10 # coordinates

A = inla.spde.make.A(mesh=mesh_sim, loc=loc)
omega_s = drop(A %*% omega_s)  #drop dimensions of an array


# Linear predictor
beta0 <- 1
sigma_e <- 0.5
lin.pred <- beta0 + omega_s

# Response variable
y <- lin.pred + sigma_e*rnorm(n)
