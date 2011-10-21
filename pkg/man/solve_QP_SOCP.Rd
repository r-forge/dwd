\name{solve_QP_SOCP}
\alias{solve_QP_SOCP}
\title{Solve a Quadratic Programming Problem}
\description{
  This routine implements the second order cone programming method from 
  Kim-Chuan Toh , Michael J. Todd, and Reha H. Tutuncu for solving quadratic 
  programming problems of the form
  \eqn{\min(-d^T b + 1/2 b^T D b)}{min(-d^T b + 1/2 b^T D b)} with the
  constraints \eqn{A^T b >= b_0}.
}
\usage{
solve_QP_SOCP(Dmat, dvec, Amat, bvec)
}
\arguments{
  \item{Dmat}{
    matrix appearing in the quadratic function to be minimized.
  }
  \item{dvec}{
    vector appearing in the quadratic function to be minimized.
  }
  \item{Amat}{
    matrix defining the constraints under which we want to minimize the
    quadratic function.
  }
  \item{bvec}{
    vector holding the values of \eqn{b_0} (defaults to zero).
  }
}
\value{
  a list with the following components:
  \item{solution}{
    vector containing the solution of the quadratic programming problem.
  }
}

\references{
Kim-Chuan Toh , Michael J. Todd, and Reha H. Tutuncu\cr
\emph{SDPT3 version 4.0 -- a MATLAB software for semidefinite-quadratic-linear 
 programming}\cr
\url{http://www.math.nus.edu.sg/~mattohkc/sdpt3.html}
}
\author{
  Hanwen Huang: \email{hanwenh@email.unc.edu};
  Perry Haaland: \email{Perry_Haaland@bd.com}; 
  Xiaosun Lu: \email{Xiaosun_Lu@bd.com}; 
  Frances Tong: \email{Frances_Tong@bd.com};
  Elaine McVey: \email{Elaine_McVey@bd.com};
  Yufeng Liu: \email{yfliu@email.unc.edu};
  J. S. Marron: \email{marron@email.unc.edu}
}

\seealso{
  \code{\link{sqlp}}
}

\examples{
##
## Assume we want to minimize: -(0 5 0) \%*\% b + 1/2 b^T b
## under the constraints:      A^T b >= b0
## with b0 = (-8,2,0)^T
## and      (-4  2  0)
##      A = (-3  1 -2)
##          ( 0  0  1)
## we can use solve.QP as follows:
##
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
solve_QP_SOCP(Dmat,dvec,Amat,bvec=bvec)
}

\keyword{optimize}
