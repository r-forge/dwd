\name{kdwd-class}
\docType{class}
\alias{kdwd-class}
\alias{scaling}
\alias{b0}
\alias{obj}
%\alias{cross,kdwd-method}
%\alias{error,kdwd-method}
%\alias{fitted,kdwd-method}
%\alias{lev,kdwd-method}
%\alias{kcall,kdwd-method}
%\alias{scaling,kdwd-method}
%\alias{type,kdwd-method}
%\alias{b0,kdwd-method}
%\alias{obj,kdwd-method}
%\alias{w,kdwd-method}
%\alias{nclass,kdwd-method}
%\alias{na.action,kdwd-method}
%\alias{index,kdwd-method}
%\alias{terms,kdwd-method}


\title{Class "kdwd" }
\description{An S4 class containing the output (model) of the
  \code{kdwd} Distance Weighted Discrimination function }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("kdwd", ...)}
  or by calls to the \code{kdwd} function.

}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Object of class \code{"character"}  containing
      the DWD type ("bdwd", "mdwd")}
    \item{\code{kcall}:}{Object of class \code{"ANY"} containing the      \code{kdwd} function call}
    \item{\code{scaling}:}{Object of class \code{"ANY"} containing the
      scaling information performed on the data}
    \item{\code{terms}:}{Object of class \code{"ANY"} containing the
      terms representation of the symbolic model used (when using a formula)}
    \item{\code{fitted}:}{Object of class \code{"output"} with the fitted values,
      predictions using the training set.}
    \item{\code{lev}:}{Object of class \code{"vector"} with the levels of the
      response (in the case of classification)}
    \item{\code{nclass}:}{Object of class \code{"numeric"}  containing
      the number of classes (in the case of classification)}
    \item{\code{w}:}{Object of class \code{"ANY"} containing the
      resulting coefficients}
    \item{\code{b0}:}{Object of class \code{"numeric"} containing the
      resulting offset }
    \item{\code{index}:}{Object of class \code{"list"} containing
      the indexes of classifiers}
      \item{\code{obj}:}{Object of class \code{vector} containing the value of the objective function. When using
      one-against-one in multiclass classification this is a vector.}
    \item{\code{error}:}{Object of class \code{"numeric"} containing the
    training error}
    \item{\code{cross}:}{Object of class \code{"numeric"} containing the
      cross-validation error }
    \item{\code{na.action}:}{Object of class \code{"ANY"} containing the
      action performed for NA }
  }
}
%\section{Methods}{
%  \describe{
%    \item{index}{\code{signature(object = "kdwd")}: return the indexes
%    of classifiers}
%    \item{cross}{\code{signature(object = "kdwd")}: returns the
%      cross-validation error }
%    \item{error}{\code{signature(object = "kdwd")}: returns the training
%      error }
%       \item{obj}{\code{signature(object = "kdwd")}: returns the value of the objective function}
%    \item{fitted}{\code{signature(object = "vm")}: returns the fitted
%      values (predict on training set) }
%    \item{lev}{\code{signature(object = "kdwd")}: returns the levels in
%      case of classification  }
%    \item{w}{\code{signature(object="kdwd")}: returns
%      the parameters of the DWD}
%    \item{kcall}{\code{signature(object="kdwd")}: returns the
%    \code{kdwd} function call}
%    \item{scaling}{\code{signature(object = "kdwd")}: returns the
%      scaling values }
%    \item{type}{\code{signature(object = "kdwd")}: returns the problem type}
%  }
%}
\author{
  Hanwen Huang: \email{hanwenh@email.unc.edu};
  Perry Haaland: \email{Perry_Haaland@bd.com}; 
  Xiaosun Lu: \email{Xiaosun_Lu@bd.com}; 
  Yufeng Liu: \email{yfliu@email.unc.edu};
  J. S. Marron: \email{marron@email.unc.edu}
}


\seealso{
  \code{\link{kdwd}}
}
\examples{
## simple example using the promotergene data set
data(promotergene)

## train a support vector machine
\donttest{gene <- kdwd(Class~.,data=promotergene,C=100,cross=4)

# the fitted values
gene@fitted
# the cross validation error
gene@cross}
}
\keyword{classes}
