% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KARMA.R
\name{KARMA}
\alias{KARMA}
\title{KARMA Model}
\usage{
KARMA(
  y,
  ar = NA,
  ma = NA,
  link = "logit",
  h = 6,
  X = NA,
  X_hat = NA,
  resid = 3
)
}
\arguments{
\item{y}{a numeric vetor or time series with UnitLindley distribution}

\item{ar}{autoregressive order}

\item{ma}{moving avarege order}

\item{link}{the link function ("logit", "probit", "cloglog")}

\item{h}{forecast horizon}

\item{X}{optionally numerical vector or matrix of external regressors, which must have the same number of rows as y.}

\item{X_hat}{the future values of \code{X} to forecast, which must have \code{h} rows.}

\item{resid}{type of residuals (1 - standardized, 2 - deviance, 3 - quartile)}
}
\value{
return a list with several informations about the model
}
\description{
Fit and forecast the Kumaraswamy Autoregressive and Moving Average (KARMA) model.
}
\examples{
x = rbeta(50, 2, 2)
fit = KARMA(ts(x), 0, 1)
fit$fitted
}
\author{
José G S Sena, Fabio M Bayer, Paulo H. Ferreira, José A Fiorucci
}
