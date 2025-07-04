% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/base.R
\name{sample.HMM}
\alias{sample.HMM}
\title{Generates an HMM process}
\usage{
sample.HMM(
  n,
  m,
  Pi = NULL,
  K = NULL,
  pp = NULL,
  emi.dist = "normal",
  emi.param = list(mu = 1:nrow(pp), sigma = rep(1, nrow(pp)))
)
}
\arguments{
\item{n}{Length of the sequence to generate}

\item{m}{Cardinality of the state space}

\item{Pi}{Initial state distribution}

\item{K}{Expected number of change points}

\item{pp}{Transition matrix}

\item{emi.dist}{Emission distribution}

\item{emi.param}{Parameters for the emission distribution}
}
\value{
An HMM process and all parameters for the estimation
}
\description{
Generates an HMM process
}
