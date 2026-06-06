if (!requireNamespace("pak", quietly = TRUE))
  install.packages("pak")
if (!requireNamespace("RTMLC", quietly = TRUE))
  pak::pkg_install("LB1304/RTMLC")
library(RTMLC)


#---- 1. Generate data ----#

# Settings
H <- 50         # Number of clusters
Enh <- 15       # Average cluster size
J <- 40         # Number of items
Jsel <- 10      # Number of items randomly selected per individual
NC.Clu <- 2     # Number of cluster-level covariates
NC.Ind <- 3     # Number of individual-level covariates
K.Clu <- 1      # Number of cluster-level latent classes
K.Ind <- 3      # Number of individual-level latent classes

# Support points
Sup.Clu.Resp <- c(0)          # Cluster-level support points for response accuracy
Sup.Ind.Resp <- c(-2, 0, 2)   # Individual-level support points for response accuracy
Sup.Clu.Time <- c(0)          # Cluster-level support points for response time
Sup.Ind.Time <- c(-1, 0, 1)   # Individual-level support points for response time

# Item parameters
dif = seq(-2, 2, length.out=J)     # Difficulties
speed0 = seq(-2, 2, length.out=J)  # Speed intercept
speed1 = seq(1, 2, length.out=J)   # Speed slope
si2 = 0.1                          # Speed variance

# Covariate parameters
Coef.Clu.Resp = seq(-1,1, length.out = NC.Clu)  # Cluster-level regression coefficients for response accuracy
Coef.Ind.Resp = seq(-1,1, length.out = NC.Ind)  # Individual-level regression coefficients for response accuracy
Coef.Clu.Time = seq(-1,1, length.out = NC.Clu)  # Cluster-level regression coefficients for response time
Coef.Ind.Time = seq(-1,1, length.out = NC.Ind)  # Individual-level regression coefficients for response time

Data <- Draw.RTMLC(
  H, Enh, J, Jsel, NC.Clu, NC.Ind,
  Sup.Clu.Resp, Sup.Ind.Resp,
  Sup.Clu.Time, Sup.Ind.Time,
  dif, speed0, speed1, si2,
  Coef.Clu.Resp, Coef.Ind.Resp,
  Coef.Clu.Time, Coef.Ind.Time
)


#---- 2. Estimate the RTM model ----#
Est0 <- Est.RTM(
  Data$Data, NC.Clu, NC.Ind, J,
  control = list(rel.tol = 1e-6)
)


#---- 3. Estimate the RTMLC model ----#
Est <- Est.RTMLC(
  Data$Data, NC.Clu, NC.Ind, J, K.Clu, K.Ind, Est0, 
  control = list(rel.tol = 1e-6)
)


#---- 4. Check estimated parameters ----#
Est$coefX.resp
Est$coefZ.resp
Est$coefX.time
Est$coefZ.time

Est$diff.resp
Est$speed0.time
Est$speed1.time
Est$si2.time

Est$sup2.resp
Est$sup2.time

table(True = Data$V, Pred = apply(Est$PostProb.Ind, 1, which.max))

