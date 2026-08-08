\name{dcsbm_fit-methods}
\alias{print.dcsbm_fit}
\alias{summary.dcsbm_fit}
\alias{print.summary.dcsbm_fit}
\alias{plot.dcsbm_fit}
\title{Methods for dcsbm_fit objects}
\description{
\code{print}, \code{summary}, and \code{plot} methods for objects of class
\code{dcsbm_fit} returned by \code{\link{fit_dcsbm}}.
}
\usage{
\method{print}{dcsbm_fit}(x, ...)

\method{summary}{dcsbm_fit}(object, ...)

\method{print}{summary.dcsbm_fit}(x, ...)

\method{plot}{dcsbm_fit}(x, type = c("loglik", "block_sizes", "heatmap"), ...)
}
\arguments{
  \item{x}{A \code{dcsbm_fit} object (or, for \code{print.summary.dcsbm_fit}, a \code{summary.dcsbm_fit} object).}
  \item{object}{A \code{dcsbm_fit} object.}
  \item{type}{For the \code{plot} method: one of \code{"loglik"} (trace plot of saved log-likelihood draws, the default), \code{"block_sizes"} (bar plot of MAP block sizes), or \code{"heatmap"} (adjacency matrix reordered by MAP block assignment).}
  \item{...}{Additional arguments; passed to the underlying base plot function for \code{plot}, otherwise unused.}
}
\value{
\code{print.dcsbm_fit} invisibly returns \code{x}.

\code{summary.dcsbm_fit} returns an object of class \code{summary.dcsbm_fit}.

\code{print.summary.dcsbm_fit} invisibly returns \code{x}.

\code{plot.dcsbm_fit} invisibly returns \code{x}; called for its plotting side effect.
}
\seealso{
\code{\link{fit_dcsbm}}
}
