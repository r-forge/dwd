\name{sqlp}
\alias{sqlp}
\title{
Solve Second Order Cone Programs
}
\description{
solve an SOCP program by infeasible path-following method
}
\usage{
sqlp(blk,At,C,b,OPTIONS,X0,y0,Z0)
}

\arguments{
  \item{blk}{a list describing the structure of the SOCP data.
}
  \item{At}{a matrix containing the coefficients for the linear and second order cone constraints.  At should have \eqn{m} columns, where \eqn{m} is the number of constraints. The number of rows in At should be \code{sum(C)}.
}
  \item{C}{a vector containing the coefficients of the objective function to be minimized.
}
  \item{b}{a vector containing the right hand side of the constraints.
}
  \item{OPTIONS}{a list that specifies parameters required in sqlp.  
}
  \item{X0}{an initial iterate for primal solution.
}
  \item{y0,Z0}{initial iterates for dual solution.
}
}
\details{
A second order cone program (SOCP) is an optimization problem similar to 
a linear program (LP), except that some variables can be constrained by 
second order cones. An exact mathematical definition can be found in 
Kim-Chuan Toh , Michael J. Todd, and Reha H. Tutuncu. This function 
implements the algorithm given in that paper which allows for multiple 
second order cone constraints as well as linear constraints. The 
objective function is given by \code{sum(C*x)} while the constraints 
are \code{A\%*\%x == b}, with \code{x} belonging to the cartesian 
product of second order cones described by \code{blk}.
}
\value{
A list containing named elements:
\item{x}{The optimal solution to the SOCP.}
\item{y,Z}{The dual solutions.}
\item{info}{Summary information.}
\item{runhist}{Run history.}
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

\examples{

#Load an example SOCP
data(sqlpData)

#Solve the socp
soln <- sqlp(blk=sqlpData$blk,At=sqlpData$At,C=sqlpData$C,b=sqlpData$b,X0=sqlpData$X0,y0=sqlpData$y0,Z0=sqlpData$Z0)

}
\keyword{optimization}
\keyword{SOCP}
