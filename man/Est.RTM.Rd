\name{Est.RTM}
\alias{Est.RTM}
\title{Estimate a response accuracy and response time model without latent classes}
\description{
Estimates a baseline item response model for the joint analysis of response
accuracy and response time without latent classes. This function can also be used
to obtain starting values for \code{\link{Est.RTMLC}}.
}
\usage{
Est.RTM(data, nc1, nc2, J, control = NULL)
}
\arguments{
  \item{data}{A data frame containing the observed data. The first two columns must be \code{cluster} and \code{individual}. They must be followed by \code{nc1} cluster-level covariates, \code{nc2} individual-level covariates, \code{J} response variables, and \code{J} response time variables. Missing responses due to the item-rotation design must be coded as \code{NA}.}
  \item{nc1}{Integer. Number of cluster-level covariates.}
  \item{nc2}{Integer. Number of individual-level covariates.}
  \item{J}{Integer. Number of items.}
  \item{control}{Optional list of control parameters. The available entries are:
  \describe{
    \item{\code{rel.tol}}{Relative tolerance for the stopping rule. Default is \code{1e-8}.}
    \item{\code{out.se}}{Logical. If \code{TRUE}, standard errors are computed. Default is \code{FALSE}.}
  }}
}
\details{
The response accuracy component is based on a normal-ogive item response model.
The response time component assumes a log-normal distribution for the response
times. The function maximizes the observed-data log-likelihood by iterating
updates for the response accuracy and response time components until convergence.

The returned estimates are arranged so that the output can be used as the
\code{param.init} argument in \code{\link{Est.RTMLC}}.
}
\value{
A list containing:
\describe{
  \item{\code{lk}}{Maximized log-likelihood.}
  \item{\code{coefX.resp}}{Estimated effects of cluster-level covariates on response accuracy, or \code{NULL} if \code{nc1 = 0}.}
  \item{\code{coefZ.resp}}{Estimated effects of individual-level covariates on response accuracy, or \code{NULL} if \code{nc2 = 0}.}
  \item{\code{diff.resp}}{Estimated item difficulty parameters.}
  \item{\code{coefX.time}}{Estimated effects of cluster-level covariates on response time, or \code{NULL} if \code{nc1 = 0}.}
  \item{\code{coefZ.time}}{Estimated effects of individual-level covariates on response time, or \code{NULL} if \code{nc2 = 0}.}
  \item{\code{speed0.time}}{Estimated item-specific time intensity parameters.}
  \item{\code{speed1.time}}{Estimated item-specific time discrimination parameters.}
  \item{\code{si2.time}}{Estimated residual variance of the log-response time distribution.}
  \item{\code{nit}}{Number of iterations.}
  \item{\code{np}}{Number of free parameters.}
  \item{\code{AIC}}{Akaike Information Criterion.}
  \item{\code{BIC}}{Bayesian Information Criterion.}
}
If \code{control$out.se = TRUE}, the list also includes standard errors for the
corresponding parameter estimates and the numerical derivative matrix
\code{Dn}.
}
\examples{
set.seed(2801)

H <- 50
Enh <- 15
J <- 20
Jsel <- 5
nc1 <- 1
nc2 <- 1

sim <- Draw.RTMLC(
  H = H, Enh = Enh, J = J, Jsel = Jsel, nc1 = nc1, nc2 = nc2,
  sup1resp = c(0), sup2resp = c(-1, 0, 1),
  sup1time = c(0), sup2time = c(0.5, 0, -0.5),
  dif = seq(-1, 1, length.out = 6),
  speed0 = rep(3, 6), speed1 = rep(1, 6), si2 = 0.2,
  coef1resp = 0.1, coef2resp = 0.2,
  coef1time = 0.1, coef2time = 0.2
)

fit0 <- Est.RTM(sim$Data, nc1 = nc1, nc2 = nc2, J = J)
fit0$lk
fit0$BIC
}
\seealso{
\code{\link{Draw.RTMLC}}, \code{\link{Est.RTMLC}}
}

