\name{kdwd}
\alias{kdwd}
\alias{kdwd,formula-method}
\alias{kdwd,matrix-method}
\title{Distance Weighted Discrimination}
\description{
  Distance Weighted Discriminations are an excellent tool for classification. 
  \code{kdwd} supports the well known binary classification formulations along with 
  native multi-class classification formulations.
}
\usage{
\S4method{kdwd}{formula}(x, data = NULL, ..., subset,na.action = na.omit, 
scaled = TRUE)
\S4method{kdwd}{matrix}(x, y = NULL, scaled = TRUE, type = "bdwd", 
C = 100, fit = TRUE, cross = 0, class.weights = NULL,
     subset, na.action = na.omit)
}

\arguments{
  \item{x}{a symbolic description of the model to be fit.  When not
    using a formula x can be a matrix containing the training data.}

  \item{data}{an optional data frame containing the training data, when using a 
  formula.}
	  
  \item{y}{a response vector with one label for each row/component of \code{x}. 
  Can be a factor.}

  \item{scaled}{A logical vector indicating the variables to be
    scaled. If \code{scaled} is of length 1, the value is recycled as
    many times as needed and all non-binary variables are scaled.}

  \item{type}{\code{kdwd} can be used for binary and multiclass classification
    Depending on whether \code{y} is a level two factor or not, the default setting 
    for \code{type} is \code{bdwd}, but can be overwritten by setting an explicit 
    value. \cr Valid options are:
    \itemize{
      \item \code{bdwd}  Binary classification or one-vs-one for multi-class 
      classification
      \item \code{mdwd}  Global multi-class classification
    }
  }

  \item{C}{cost of constraints violation (default: 1) this is the
    `C'-constant of the regularization term in the Lagrange
    formulation.}

  \item{class.weights}{a named vector of weights for the different
    classes, used for asymmetric class sizes. Not all factor levels have
    to be supplied (default weight: 1). All components have to be named.}

  \item{cross}{if a integer value k>0 is specified, a k-fold cross
    validation on the training data is performed to assess the quality
    of the model: the accuracy rate for classification.}

  \item{fit}{indicates whether the fitted values should be computed
    and included in the model or not (default: \code{TRUE}).}

  \item{subset}{An index vector specifying the cases to be used in the
          training sample.  (NOTE: If given, this argument must be
          named.)}

   \item{na.action}{A function to specify the action to be taken if \code{NA}s are
          found. The default action is \code{na.omit}, which leads to rejection of 
          cases with missing values on any required variable. An alternative
	  is \code{na.fail}, which causes an error if \code{NA} cases
	  are found. (NOTE: If given, this argument must be named.)}
  \item{\dots}{further arguments.}
        }
    
    \value{
      An S4 object of class \code{"kdwd"} containing the fitted model, Accessor 
      functions can be used to access 
      the slots of the object (see examples) which include:
  \item{w}{The direction vector pointing towards positive class, unit vector 
  (i.e length 1). In case of	one-against-one classification this is a list of vectors. 
  In case of global multiclass classification,
  	it is a matrix.}
  \item{b0}{The intercept. In case of one-against-one classification this is a 
  list of scalers. In case of 
  	global multiclass classification, it is a vector.}
  \item{obj}{The value of the objective function. In case of one-against-one 
  classification this is a vector 
  	of values.} 
  \item{error}{Training error.}
  \item{cross}{Cross validation error, (when cross > 0).}
  \item{fitted}{Fitted class label for training data.}
  }

  
\details{
  \code{kdwd} uses SDPT3 infeasible path-following algorithm to solve the DWD SOCP 
  problem. \cr
  For multiclass-classification with \eqn{k} classes, \eqn{k > 2}, \code{bdwd} uses 
  the
  `one-against-one'-approach, in which \eqn{k(k-1)/2} binary classifiers are
  trained; the appropriate class is found by a voting scheme,
  The \code{mdwd} formulations deal with the multiclass-classification problems by 
  solving 
  a single SOCP problem involving all the classes.\cr
  In classification when \code{cross} is \code{k} a k-fold cross validation is
  performed on the data.
  The data can be passed to the \code{kdwd} function in a \code{matrix} or a
  \code{data.frame}.
}
\references{
  \itemize{
    \item
      Kim-Chuan Toh , Michael J. Todd, and Reha H. Tutuncu\cr
      \emph{SDPT3 version 4.0 -- a MATLAB software for semidefinite-quadratic-linear 
      programming}\cr
      \url{http://www.math.nus.edu.sg/~mattohkc/sdpt3.html}
   \item
      J. S. Marron and Michael Todd\cr
       \emph{Distance Weighted Discrimination}
       \url{http://ecommons.cornell.edu/bitstream/1813/9217/1/TR001339.pdf}
  }
}
\author{
  Hanwen Huang: \email{hanwenh@email.unc.edu};
  Perry Haaland: \email{Perry_Haaland@bd.com}; 
  Xiaosun Lu: \email{Xiaosun_Lu@bd.com}; 
  Yufeng Liu: \email{yfliu@email.unc.edu};
  J. S. Marron: \email{marron@email.unc.edu}
}
\seealso{\code{\link{predict.kdwd}}, \code{\link{kdwd-class}} }
  
\keyword{methods}
\keyword{classification}
\keyword{optimization}

\examples{
    ## simple binary classification example
	data(promotergene)

	## train DWD
	gene <- kdwd(Class~.,data=promotergene,C=100,scaled=TRUE,cross=5)

	gene@fitted

	## simple multiclass example using the famous iris data
	data(iris)

	## train an OVO multiclass DWD
	\donttest{irismodel <- kdwd(Species~.,data=iris,type="bdwd",C=100,scaled=TRUE,cross=5)

	## get fitted values
	irismodel@fitted

	## Test on the training set
	predict(irismodel, iris)}
}
