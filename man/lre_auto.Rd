% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lrem.R
\name{lre_auto}
\alias{lre_auto}
\title{A combination of the Klein and Blanchard, Kahn methods.}
\usage{
lre_auto(A, E = NULL, n)
}
\arguments{
\item{A}{Coefficient matrix on x(t)}

\item{E}{Coefficient matrix on x(t+1), if not given, the function evaluates the BK method}

\item{n}{number of predetermined variables}
}
\value{
Motion and decision rule functions according to either the Blanchard Kahn approach or through the Klein method
}
\description{
A combination of the Klein and Blanchard, Kahn methods.
}
