% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_prob.R
\name{argH3}
\alias{argH3}
\title{Compute best path with three segments and jump at k1,k2 on the interval l:r}
\usage{
argH3(l, k1, k2, r, x0, par)
}
\arguments{
\item{l}{Left endpoint of the interval}

\item{k1}{First change point}

\item{k2}{Second change point}

\item{r}{Right endpoint of the interval}

\item{x0}{Previous state}

\item{par}{Model parameters}
}
\value{
Change points (\code{k1, k2}), state values and H3-value
}
\description{
Compute best path with three segments and jump at k1,k2 on the interval l:r
}
