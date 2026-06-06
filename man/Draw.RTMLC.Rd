\name{Draw.RTMLC}
\alias{Draw.RTMLC}
\title{Simulate data from a Response Time Multilevel Latent Class IRT model}
\description{
Generates clustered item-response data from a Response Time Multilevel Latent
Class (RTMLC) IRT model. The generated data include dichotomous item responses,
item response times, cluster-level covariates, individual-level covariates, and
the true latent class allocations.
}
\usage{
Draw.RTMLC(H, Enh, J, Jsel, nc1, nc2,
           sup1resp, sup2resp, sup1time, sup2time,
           dif, speed0, speed1, si2,
           coef1resp, coef2resp, coef1time, coef2time)
}
\arguments{
  \item{H}{Integer. Number of clusters.}
  \item{Enh}{Numeric. Expected cluster size used as the mean of the Poisson distribution from which cluster sizes are generated. Cluster sizes are constrained to be at least one.}
  \item{J}{Integer. Number of items in the item pool.}
  \item{Jsel}{Integer. Number of items randomly selected and administered to each individual.}
  \item{nc1}{Integer. Number of cluster-level covariates.}
  \item{nc2}{Integer. Number of individual-level covariates.}
  \item{sup1resp}{Numeric vector. Cluster-level support points for the response accuracy component. Its length defines the number of cluster-level latent classes.}
  \item{sup2resp}{Numeric vector. Individual-level support points for the response accuracy component. Its length defines the number of individual-level latent classes.}
  \item{sup1time}{Numeric vector. Cluster-level support points for the response time component.}
  \item{sup2time}{Numeric vector. Individual-level support points for the response time component.}
  \item{dif}{Numeric vector of length \code{J}. Item difficulty parameters for the response accuracy component.}
  \item{speed0}{Numeric vector of length \code{J}. Item-specific time intensity parameters.}
  \item{speed1}{Numeric vector of length \code{J}. Item-specific time discrimination parameters.}
  \item{si2}{Numeric. Residual variance of the log-response time distribution.}
  \item{coef1resp}{Numeric vector of length \code{nc1}. Effects of cluster-level covariates on response accuracy.}
  \item{coef2resp}{Numeric vector of length \code{nc2}. Effects of individual-level covariates on response accuracy.}
  \item{coef1time}{Numeric vector of length \code{nc1}. Effects of cluster-level covariates on response time.}
  \item{coef2time}{Numeric vector of length \code{nc2}. Effects of individual-level covariates on response time.}
}
\details{
For each cluster, a cluster-level latent class is sampled uniformly from the
available cluster-level classes. For each individual, an individual-level latent
class is sampled uniformly from the available individual-level classes. Cluster-
and individual-level covariates are generated from independent standard normal
distributions.

For each individual, \code{Jsel} items are selected at random from the item pool.
Responses are generated from a normal-ogive response accuracy model, while
response times are generated from a log-normal response time model. Items not
administered to a given individual are coded as \code{NA} in both the response and
time matrices.
}
\value{
A list with the following components:
\describe{
  \item{\code{Data}}{A data frame containing the simulated data. The columns are
  \code{cluster}, \code{individual}, the cluster-level covariates \code{X1}, ...,
  the individual-level covariates \code{Z1}, ..., the response variables
  \code{resp1}, ..., \code{respJ}, and the response time variables
  \code{time1}, ..., \code{timeJ}.}
  \item{\code{U}}{Integer vector of length \code{H} containing the true cluster-level latent class allocations.}
  \item{\code{V}}{Integer vector of length equal to the total number of individuals containing the true individual-level latent class allocations.}
}
}
\examples{
set.seed(2801)

H <- 50
Enh <- 15
J <- 20
Jsel <- 5
nc1 <- 1
nc2 <- 1

out <- Draw.RTMLC(
  H = H, Enh = Enh, J = J, Jsel = Jsel, nc1 = nc1, nc2 = nc2,
  sup1resp = c(0), sup2resp = c(-1, 0, 1),
  sup1time = c(0), sup2time = c(0.5, 0, -0.5),
  dif = seq(-1, 1, length.out = J),
  speed0 = rep(3, J), speed1 = rep(1, J), si2 = 0.2,
  coef1resp = 0.1, coef2resp = 0.2,
  coef1time = 0.1, coef2time = 0.2
)

head(out$Data)
table(out$V)
}
\seealso{
\code{\link{Est.RTM}}, \code{\link{Est.RTMLC}}
}

