## Distance Weight Discrimination
## author : Hanwen Huang
## updated : 08.04.11

setClass("kdwd", representation(scaling = "ANY",
				type = "character",
				fitted = "ANY",  
				lev = "vector",
				nclass = "numeric",
				error = "vector",
				cross = "vector",
				na.action= "ANY",
				terms = "ANY",
				w = "ANY",
				b0 = "ANY",
				obj = "ANY",
				index = "ANY",
				kcall = "call"
		))

if(!isGeneric("kcall")){
  if (is.function("kcall"))
    fun <- kcall
  else fun <- function(object) standardGeneric("kcall")
  setGeneric("kcall", fun)
}
setMethod("kcall", "kdwd", function(object) object@kcall)
setGeneric("kcall<-", function(x, value) standardGeneric("kcall<-"))
setReplaceMethod("kcall", "kdwd", function(x, value) {
  x@kcall <- value
  x
})

if(!isGeneric("type")){
  if (is.function("type"))
    fun <- type
  else fun <- function(object) standardGeneric("type")
  setGeneric("type", fun)
}
setMethod("type", "kdwd", function(object) object@type)
setGeneric("type<-", function(x, value) standardGeneric("type<-"))
setReplaceMethod("type", "kdwd", function(x, value) {
  x@type <- value
  x
})

setMethod("terms", "kdwd", function(x, ...) x@terms)
setGeneric("terms<-", function(x, value) standardGeneric("terms<-"))
setReplaceMethod("terms", "kdwd", function(x, value) {
  x@terms <- value
  x
})

setMethod("fitted", "kdwd", function(object, ...) object@fitted)
setGeneric("fitted<-", function(x, value) standardGeneric("fitted<-"))
setReplaceMethod("fitted", "kdwd", function(x, value) {
  x@fitted <- value
  x
})

if(!isGeneric("lev")){
  if (is.function("lev"))
    fun <- lev
  else fun <- function(object) standardGeneric("lev")
  setGeneric("lev", fun)
}
setMethod("lev", "kdwd", function(object) object@lev)
setGeneric("lev<-", function(x, value) standardGeneric("lev<-"))
setReplaceMethod("lev", "kdwd", function(x, value) {
  x@lev <- value
  x
})

if(!isGeneric("nclass")){
  if (is.function("nclass"))
    fun <- nclass
  else fun <- function(object) standardGeneric("nclass")
  setGeneric("nclass", fun)
}
setMethod("nclass", "kdwd", function(object) object@nclass)
setGeneric("nclass<-", function(x, value) standardGeneric("nclass<-"))
setReplaceMethod("nclass", "kdwd", function(x, value) {
  x@nclass <- value
  x
})

if(!isGeneric("index")){
  if (is.function("index"))
    fun <- index
  else fun <- function(object) standardGeneric("index")
  setGeneric("index", fun)
}
setMethod("index", "kdwd", function(object) object@index)
setGeneric("index<-", function(x, value) standardGeneric("index<-"))
setReplaceMethod("index", "kdwd", function(x, value) {
  x@index <- value
  x
})

if(!isGeneric("error")){
  if (is.function("error"))
    fun <- error
  else fun <- function(object) standardGeneric("error")
  setGeneric("error", fun)
}
setMethod("error", "kdwd", function(object) object@error)
setGeneric("error<-", function(x, value) standardGeneric("error<-"))
setReplaceMethod("error", "kdwd", function(x, value) {
  x@error <- value
  x
})

if(!isGeneric("cross")){
  if (is.function("cross"))
    fun <- cross
  else fun <- function(object) standardGeneric("cross")
  setGeneric("cross", fun)
}
setMethod("cross", "kdwd", function(object) object@cross)
setGeneric("cross<-", function(x, value) standardGeneric("cross<-"))
setReplaceMethod("cross", "kdwd", function(x, value) {
  x@cross <- value
  x
})

if(!isGeneric("na.action")){
  if (is.function("na.action"))
    fun <- na.action
  else fun <- function(object) standardGeneric("na.action")
  setGeneric("na.action", fun)
}
setMethod("na.action", "kdwd", function(object) object@na.action)
setGeneric("na.action<-", function(x, value) standardGeneric("na.action<-"))
setReplaceMethod("na.action", "kdwd", function(x, value) {
  x@na.action <- value
  x
})

if(!isGeneric("scaling")){
  if (is.function("scaling"))
    fun <- scaling
  else fun <- function(object) standardGeneric("scaling")
  setGeneric("scaling", fun)
}
setMethod("scaling", "kdwd", function(object) object@scaling)
setGeneric("scaling<-", function(x, value) standardGeneric("scaling<-"))
setReplaceMethod("scaling", "kdwd", function(x, value) {
  x@scaling<- value
  x
})

if(!isGeneric("obj")){
  if (is.function("obj"))
    fun <- obj
  else fun <- function(object) standardGeneric("obj")
  setGeneric("obj", fun)
}
setMethod("obj", "kdwd", function(object) object@obj)
setGeneric("obj<-", function(x, value) standardGeneric("obj<-"))
setReplaceMethod("obj", "kdwd", function(x, value) {
  x@obj<- value
  x
})


if(!isGeneric("w")){
  if (is.function("w"))
    fun <- w
  else fun <- function(object) standardGeneric("w")
  setGeneric("w", fun)
}
setMethod("w", "kdwd", function(object) object@w)
setGeneric("w<-", function(x, value) standardGeneric("w<-"))
setReplaceMethod("w", "kdwd", function(x, value) {
  x@w <- value
  x
})

if(!isGeneric("b0")){
  if (is.function("b0"))
    fun <- b0
  else fun <- function(object) standardGeneric("b0")
  setGeneric("b0", fun)
}
setMethod("b0", "kdwd", function(object) object@b0)
setGeneric("b0<-", function(x, value) standardGeneric("b0<-"))
setReplaceMethod("b0", "kdwd", function(x, value) {
  x@b0 <- value
  x
})

setGeneric("kdwd", function(x, ...) standardGeneric("kdwd"))
setMethod("kdwd",signature(x="formula"),
          function (x, data=NULL, ..., subset, na.action = na.omit,
                    scaled = TRUE){
            cl <- match.call()
            m <- match.call(expand.dots = FALSE)
            if (is.matrix(eval(m$data, parent.frame())))
              m$data <- as.data.frame(data)
            m$... <- NULL
            m$formula <- m$x
            m$x <- NULL
            m$scaled <- NULL
            m[[1]] <- as.name("model.frame")
            m <- eval(m, parent.frame())
            Terms <- attr(m, "terms")
            attr(Terms, "intercept") <- 0    ## no intercept
            x <- model.matrix(Terms, m)
            y <- model.extract(m, response)
            if (length(scaled) == 1)
              scaled <- rep(scaled, ncol(x))
            if (any(scaled)) {
              remove <- unique(c(which(labels(Terms) %in% names(attr(x, "contrasts"))),
                                 which(!scaled)
                                 )
                               )
              scaled <- !attr(x, "assign") %in% remove
            }
            ret <- kdwd(x, y, scaled = scaled, ...)
            kcall(ret) <- cl
            attr(Terms,"intercept") <- 0 ## no intercept
            terms(ret) <- Terms
            if (!is.null(attr(m, "na.action")))
              na.action(ret) <- attr(m, "na.action")
            return (ret)
          })

setMethod("kdwd",signature(x="matrix"),
          function (x,
                    y         = NULL,
                    scaled    = TRUE,
                    type      = "bdwd",
                    C         = 100,
                    fit       = TRUE,
                    cross     = 0,
                    class.weights = NULL,
                    subset, 
                    na.action = na.omit)
          { 
            ret <- new("kdwd")
            if (!missing(subset)) 
              x <- x[subset,]
            if (is.null(y))
              x <- na.action(x)
            else {
              df <- na.action(data.frame(y, x))
              y <- df[,1]
              x <- as.matrix(df[,-1])
            }
            ncols <- ncol(x)
            m <- nrows <- nrow(x)
            na.action(ret) <- na.action
            
            if(!is.null(type))
              type(ret) <- match.arg(type,c("bdwd","mdwd"))
            
            x.scale <- NULL
            ## scaling
            if (length(scaled) == 1){
              scaled <- rep(scaled, ncol(x))
            }
            if (any(scaled)) {
              co <- !apply(x[,scaled, drop = FALSE], 2, var)
              if (any(co)) {
                scaled <- rep(FALSE, ncol(x))
                warning(paste("Variable(s)",
                              paste("`",colnames(x[,scaled, drop = FALSE])[co],
                                    "'", sep="", collapse=" and "),
                              "constant. Cannot scale data.")
                        )
              } else {
                xtmp <- scale(x[,scaled])
                x[,scaled] <- xtmp
                x.scale <- attributes(xtmp)[c("scaled:center","scaled:scale")]
              }
            }
            
            if (!is(y,"vector") && !is.factor (y)) 
              stop("y must be a vector or a factor.")
            
            weightlabels <- NULL
            weight <- NULL
            ## in case of classification: transform factors into integers
            if (is.factor(y)) {
              lev(ret) <- levels(y)
              ybuf = y		
              y <- as.integer(y)
            }else {
              lev(ret) <- sort(unique(y))
            }
            if (!is.null(class.weights)) {
              weightlabels <- match(names(class.weights),lev(ret))
              if (any(is.na(weightlabels)))
                stop ("At least one level name is missing or misspelled.")
            }
            ## initialize    
            nclass(ret) <- length(unique(y))
            p <- 0
            K <- 0 
            problem <- NULL
            
            ## C classification
            if(type(ret) == "bdwd"){
              indexes <- lapply(sort(unique(y)), function(kk) which(y == kk))
              for (i in 1:(nclass(ret)-1)) {
                jj <- i+1
                for(j in jj:nclass(ret)) {
                  p <- p+1
                  ## prepare the data
                  if(!is.null(class.weights)){
                    weight <- class.weights[weightlabels[c(i,j)]]
                  }
                  trainp = t(x[indexes[[i]],])
                  trainn = t(x[indexes[[j]],])
                  dwd.fit <- DWD1SM_weight(trainp,trainn,threshfact=C,
                                           wgt=weight)
                  
                  w(ret)[p] <- list(dwd.fit$w)
                  b0(ret)[p] <- list(dwd.fit$beta)
                  obj(ret)[p] <- list(dwd.fit$obj)
                  index(ret)[p] <- list(c(i,j))
                }
              }
            } 
            
            ## global multiclass classification  
            if(type(ret) =="mdwd"){
              if(!is.null(class.weights))
                weight <- class.weights[weightlabels]
              dwd.fit <- mdwd(t(x),y,threshfact=C,wgt=weight)
              w(ret) <- dwd.fit$w
              b0(ret) <- dwd.fit$beta
              obj(ret) <-  dwd.fit$obj
              index(ret) <- NULL
            }
            kcall(ret) <- match.call()
            
            fitted(ret)  <- if (fit) predict(ret, x) else NULL
            
            
            if(any(scaled))
              scaling(ret) <- list(scaled = scaled, x.scale = x.scale)
            else
              scaling(ret) <- NULL
            
            if (fit){
              error(ret) <- 1 - .classAgreement(table(as.factor(ybuf),fitted(ret)))
            }
            
            cross(ret) <- -1
            if(cross == 1)
              cat("\n","cross should be >1 no cross-validation done!","\n","\n")
            else if (cross > 1)
              {
                cerror <- 0
                suppressWarnings(vgr<-split(sample(1:m,m),1:cross))
                for(i in 1:cross)
                  {
                    cind <-  unsplit(vgr[-i],factor(rep((1:cross)[-i],unlist(lapply(vgr[-i],length)))))
                    if(is.null(class.weights))
                      cret <- kdwd(x[cind,],y[cind],type=type(ret),C=C,scaled=FALSE,
                                   cross=0,fit=FALSE)
                    else
                      cret <- kdwd(x[cind,],y[cind],type=type(ret),C=C,scaled=FALSE, 
                                   cross=0,fit=FALSE,class.weights=class.weights)
                    cres <- predict(cret, x[vgr[[i]],,drop=FALSE])
                    cerror <- (1 - .classAgreement(table(y[vgr[[i]]],as.integer(cres))))/cross + cerror
                  }
                cross(ret) <- cerror
              }
            return(ret)
          })


.classAgreement <- function (tab) {
  n <- sum(tab)
  if (!is.null(dimnames(tab))) {
    lev <- intersect(colnames(tab), rownames(tab))
    p0 <- sum(diag(tab[lev, lev])) / n
  } else {
    m <- min(dim(tab))
    p0 <- sum(diag(tab[1:m, 1:m])) / n
  }
  return(p0)
}

##**************************************************************#
## predict for matrix, data.frame input

setMethod("predict", signature(object = "kdwd"),
          function (object, newdata)
          {
            
            if(!is(newdata,"list")){ 
              if (!is.null(terms(object)))
                {
                  if(!is.matrix(newdata))
                    newdata <- model.matrix(delete.response(terms(object)), 
                                            as.data.frame(newdata), na.action = na.action(object))
                }
              else
                newdata  <- if (is.vector(newdata)) t(t(newdata)) else as.matrix(newdata)
              
              
              newnrows <- nrow(newdata)
              newncols <- ncol(newdata)
            }
            else
              newnrows <- length(newdata)
            
            p <- 0
            
            if (is.list(scaling(object))){				
              newdata[,scaling(object)$scaled] <-
                scale(newdata[,scaling(object)$scaled, drop = FALSE],
                      center = scaling(object)$x.scale$"scaled:center", 
                      scale  = scaling(object)$x.scale$"scaled:scale")
            }
            if(type(object) == "bdwd")
              {
                predres <- 1:newnrows
                votematrix <- matrix(0,nclass(object),newnrows)
                
                for(i in 1:(nclass(object)-1))
                  {
                    jj <- i+1
                    for(j in jj:nclass(object))
                      {
                        p <- p+1
                        
                        ret <- newdata%*%w(object)[[p]] + b0(object)[[p]][1,1]
                        votematrix[i,ret>0] <- votematrix[i,ret>0] + 1
                        votematrix[j,ret<0] <- votematrix[j,ret<0] + 1
                      }
                  }
                predres <- sapply(predres, function(x) which.max(votematrix[,x]))
              }
            
            if(type(object) == "mdwd")
              {
                predres <- 1:newnrows
                votematrix <- matrix(0,nclass(object),newnrows)
                for(i in 1:nclass(object)){
                  votematrix[i,] <- newdata%*%w(object)[,i]+b0(object)[i]
                }
                predres <- sapply(predres, function(x) which.max(votematrix[,x]))
              }
            predClass = lev(object)[predres]
            return(predClass)
          })

