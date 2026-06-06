\name{Est.RTMLC}
\alias{Est.RTMLC}
\title{Estimate a Response Time Multilevel Latent Class IRT model}
\description{
Estimates a Response Time Multilevel Latent Class (RTMLC) IRT model for the
joint analysis of response accuracy and response time in multilevel item-response
data.
}
\usage{
Est.RTMLC(data, nc1, nc2, J, k1, k2, param.init, control = NULL)
}
\arguments{
  \item{data}{A data frame containing the observed data. The first two columns must be \code{cluster} and \code{individual}. They must be followed by \code{nc1} cluster-level covariates, \code{nc2} individual-level covariates, \code{J} response variables, and \code{J} response time variables. Missing responses due to the item-rotation design must be coded as \code{NA}.}
  \item{nc1}{Integer. Number of cluster-level covariates.}
  \item{nc2}{Integer. Number of individual-level covariates.}
  \item{J}{Integer. Number of items.}
  \item{k1}{Integer. Number of cluster-level latent classes.}
  \item{k2}{Integer. Number of individual-level latent classes.}
  \item{param.init}{List of initial parameter values. It is typically obtained from \code{\link{Est.RTM}}. The list must contain the entries \code{coefX.resp}, \code{coefZ.resp}, \code{diff.resp}, \code{speed0.time}, \code{speed1.time}, \code{coefX.time}, \code{coefZ.time}, and \code{si2.time}.}
  \item{control}{Optional list of control parameters. The available entries are:
  \describe{
    \item{\code{rel.tol}}{Relative tolerance for the stopping rule. Default is \code{1e-8}.}
    \item{\code{out.se}}{Logical. If \code{TRUE}, standard errors are computed. Default is \code{FALSE}.}
    \item{\code{init.rand}}{Logical. If \code{TRUE}, random starting values are used for the support points and class weights. Default is \code{FALSE}.}
  }}
}
\details{
The model includes discrete latent variables at both the cluster and individual
levels. Each latent class is characterized by support points for response
accuracy and response time. The response accuracy component uses a normal-ogive
IRT parameterization, whereas the response time component assumes a conditional
log-normal distribution.

The function implements maximum likelihood estimation through an
Expectation-Maximization algorithm. The E-step computes posterior class
probabilities for the cluster- and individual-level latent variables. The M-step
updates the class weights, the response accuracy parameters, the response time
parameters, and the residual variance.

The first support point at each level is fixed to zero for identification.
Similarly, the first time discrimination parameter is used as the reference
parameter in the response time component.
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
  \item{\code{piv1}}{Estimated weights of the cluster-level latent classes.}
  \item{\code{piv2}}{Estimated weights of the individual-level latent classes.}
  \item{\code{sup1.resp}}{Estimated cluster-level support points for response accuracy.}
  \item{\code{sup2.resp}}{Estimated individual-level support points for response accuracy.}
  \item{\code{sup1.time}}{Estimated cluster-level support points for response time.}
  \item{\code{sup2.time}}{Estimated individual-level support points for response time.}
  \item{\code{nit}}{Number of EM iterations.}
  \item{\code{AIC}}{Akaike Information Criterion.}
  \item{\code{BIC}}{Bayesian Information Criterion.}
}
If \code{control$out.se = TRUE}, the list also includes standard errors for the
corresponding parameter estimates, standard errors for the transformed latent
class weight parameters, and the numerical derivative matrix \code{Dn}.
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

init <- Est.RTM(sim$Data, nc1 = nc1, nc2 = nc2, J = J)

fit <- Est.RTMLC(
  data = sim$Data, nc1 = nc1, nc2 = nc2, J = J,
  k1 = 1, k2 = 3, param.init = init
)

fit$lk
fit$piv2
fit$sup2.resp
fit$sup2.time
}
\seealso{
\code{\link{Draw.RTMLC}}, \code{\link{Est.RTM}}
}

