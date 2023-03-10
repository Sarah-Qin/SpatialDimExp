% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spade.R
\name{GenerateData}
\alias{GenerateData}
\title{Generate three dimensional geographical data.}
\usage{
GenerateData(n, SimRep, trueTheta)
}
\arguments{
\item{n}{The number of locations.}

\item{SimRep}{The number of realizations of the Gaussian process \eqn{Y(X)}.}

\item{trueTheta}{A given positive value that controls the correlation of the spatial random process.}
}
\value{
\itemize{
 \item LocIndex: three dimensional location coordinates;
 \item YObs: realizations of the Gaussian process \eqn{Y(X)};
 \item SigmaObs: covariance matrix of the random process \eqn{Y(X)}.
}
}
\description{
The spatial random process \eqn{Y(X)} follows a Gaussian process.
The locations coordinates are simulated on a three-dimensional half-ellipsoid
centered at (0,0,0) and the projection of the first two dimensions
is a disk centered at the origin.
}
\references{
\itemize{
 \item Bornn, L.,Shaddick, G., and Zidek, J.V.(2012).
   Modeling non-stationary processes through dimension expansion.
   Journal of American Statistical Association, 107:281–289.
 \item Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020).
   Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
}
}
