\name{sqlpData}
\alias{sqlpData}
\title{SOCP input Data}
\description{A simulated example which can be directed input into \code{sqlp}.}
\usage{data(sqlpData)}

\format{
  A list of 7 entries which serve as the direct input arguments for \code{sqlp}.
}
\details{
  This dataset contains a list of 7 entries: 
  \itemize{
  \item{blk}{a list describing the structure of the SOCP data.}
  \item{At}{a matrix containing the coefficients for the linear and second order cone constraints.  At should have \eqn{m} columns, where \eqn{m} is the number of constraints. The number of rows in At should be \code{sum(C)}.}
  \item{C}{a vector containing the coefficients of the objective function to be minimized.}
  \item{b}{a vector containing the right hand side of the constraints.}
  \item{X0}{an initial iterate for primal solution.}
  \item{y0,Z0}{initial iterates for dual solution.}
  }
  }
\keyword{datasets}