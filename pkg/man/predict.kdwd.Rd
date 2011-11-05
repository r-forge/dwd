\name{predict.kdwd}
\alias{predict.kdwd}
\alias{predict,kdwd-method}
\title{predict method for DWD object}

\description{Prediction of test data using distance weighted discrimination}

\usage{
\S4method{predict}{kdwd}(object, newdata)
}

\arguments{

  \item{object}{an S4 object of class \code{kdwd} created by the
    \code{kdwd} function}
  \item{newdata}{a data frame or matrix containing new data}
}

\value{
predicted classes.
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

\keyword{methods}
\keyword{classification}
\keyword{optimization}

\examples{

## example using the promotergene data set
data(promotergene)

## create test and training set
ind <- sample(1:dim(promotergene)[1],20)
genetrain <- promotergene[-ind, ]
genetest <- promotergene[ind, ]

## train a distance weighted discrimination
\donttest{gene <- kdwd(Class~.,data=genetrain,C=100,cross=5,scaled=TRUE)
gene@fitted

## predict gene type on the test set
genetype <- predict(gene,genetest)
genetype}
}
