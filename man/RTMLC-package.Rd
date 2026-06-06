\name{RTMLC-package}
\alias{RTMLC-package}
\alias{RTMLC}
\docType{package}

\title{Response Time Multilevel Latent Class IRT Models}

\description{
The \pkg{RTMLC} package provides functions to simulate and estimate Response
Time Multilevel Latent Class (RTMLC) IRT models for the joint analysis of
response accuracy and response time in multilevel item-response data.
}

\details{
The package is designed for dichotomously scored items administered to
individuals nested within higher-level units, such as students within school
classes. The RTMLC model relies on discrete latent variables at both the
individual and cluster levels, allowing unobserved heterogeneity to be
represented through latent classes.

Response accuracy is modeled through a normal-ogive IRT sub-model, whereas
response times are modeled by a log-normal sub-model. Individual- and
cluster-level covariates can be included in both components of the model. The
estimation procedure is based on maximum likelihood and is implemented through
an Expectation-Maximization algorithm.

The main user-level functions are:
\itemize{
  \item \code{\link{Draw.RTMLC}}, to simulate clustered item-response data with
  response accuracy and response times under an RTMLC model;
  \item \code{\link{Est.RTM}}, to estimate a baseline response accuracy and
  response time model without latent classes, which can also be used to obtain
  starting values;
  \item \code{\link{Est.RTMLC}}, to estimate the full Response Time Multilevel
  Latent Class IRT model.
}

The data frame required by the estimation functions must contain, in this order:
cluster identifiers, individual identifiers, cluster-level covariates,
individual-level covariates, item responses, and item response times. Missing
responses due to the item-rotation design must be coded as \code{NA}.
}

\author{
Luca Brusa \email{luca.brusa@unimib.it},
Francesco Bartolucci,
Fulvia Pennoni,
Giorgio Vittadini.

Maintainer: Luca Brusa \email{luca.brusa@unimib.it}
}

\references{
Brusa, L., Bartolucci, F., Pennoni, F., and Vittadini, G. Joint modeling of
response accuracy and time in process data: A multilevel latent class IRT model.
}

\seealso{
\code{\link{Draw.RTMLC}},
\code{\link{Est.RTM}},
\code{\link{Est.RTMLC}}
}

\examples{
## See the help pages of the main functions for examples:
help("Draw.RTMLC")
help("Est.RTM")
help("Est.RTMLC")
}

\keyword{package}
\keyword{models}
