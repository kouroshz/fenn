#' Free Energy Nearest Nighbour (FENN)
#'
#' @description Fits a Von Neumann entropy penalized distance metric learning model.
#'
#' @usage fenn(x, \dots)
#' 
#' \method{fenn}{formula}(formula, data, \dots, subset, na.action)
#'
#' \method{fenn}{default}(x, grouping, prior = proportions, tol = 1.0e-4, \dots)
#'                       
#' \method{fenn}{data.frame}(x, \dots)
#'
#' \method{fenn}{matrix}(x, grouping, \dots, subset, na.action)
#'
#' @param formula A formula of the form \code{groups ~ x1 + x2 + \dots}  That is, the
#'                response is the grouping factor and the right hand side specifies
#'                the (non-factor) discriminators.
#'
#' @param data Data frame from which variables specified in \code{formula} are
#'             preferentially to be taken.
#'
#' @param x (required if no formula is given as the principal argument.)
#'           a matrix or data frame or Matrix containing the explanatory variables.
#'
#' @param grouping (required if no formula principal argument is given.)
#'                 a factor specifying the class for each observation.
#'
#' @param prior the prior probabilities of class membership.  If unspecified, the
#'              class proportions for the training set are used.  If present, the
#'              probabilities should be specified in the order of the factor levels.
#'
#' @param tol A tolerance to decide if a matrix is singular; it will be used to modify
#'            the scatter matrices by stabilizing eignevalues less that \code{tol^2}.
#'
#' @param subset An index vector specifying the cases to be used in the training
#'               sample.  (NOTE: If given, this argument must be named.)
#'
#' @param na.action A function to specify the action to be taken if \code{NA}s are found.
#'                  The default action is for the procedure to fail.  An alternative is
#'                  \code{na.omit}, which leads to rejection of cases with missing values on
#'                  any required variable.  (NOTE: If given, this argument must be named.)
#'
#' @param \dots arguments passed to or from other methods.
#'
#' @return An object of class \code{"fenn"} containing the following components:
#' \item{prior}{The prior probabilities used.}
#' \item{counts}{The group counts.}
#' \item{means}{The group means.}
#' \item{X.S}{The matrix of projections into similarity directions. Same as average within class
#'            covariance matrices.}
#' \item{X.D}{The matrix of projections into dissimilarity directions. Same as \code{X.S + (n/n-1) * Cov(m)},
#'            where \code{m} is the class means.}
#' \item{X_D_neg_1_2}{The matrix \code{X.D} raised to the power -1/2. This is the scaling that is applied
#'                    to data points prior to tilde transformation.}
#' \item{S_1_2}{The learned optimal transformation that should be applied to tilde transformed data points.}
#' \item{X_tilde_S}{The hamiltonian generated from data points.}
#' \item{informative.dims}{The index of informative directions in the optimal space. Can be used for dimension reduction.}
#' \item{muVec}{The automatically selected path of \code{mu} (temperature) values.}
#' \item{best.mu}{The the optimal (temperature) parameter obtained by maximizing Fisher Information.}
#' \item{E}{The average Energy for an automatically selected path of \code{mu} values \code{muVec}.}
#' \item{dE}{The Fisher Information of \code{mu}.}
#' \item{scaling}{The weights (eigenvalues) of maximum dissimilarity directions (eigenvectors of \code{S_1_2}).}
#' \item{x.tilde}{The tilde transformed data.}
#' \item{x.fenn}{The \code{fenn} transformed data.}
#' \item{N}{The number of observations used.}
#' \item{groupings}{The class variable of original data points.}
#' \item{call}{The (matched) function call.}
#' 
#' @details The function fits a Von Neumann Entropy penalized distance metric learning problem
#' to identify informative features and directions of maximam dissimilarity in the multi class case.
#' The method automatically finds the optimal value of the entropy tuning parameter by maximizing the
#' Fisher Information. The method can be used for sinlge class case to identify informative directions 
#' as well as multi class case to identiy directions of maximum dissimilarity. In the multi class case,
#' optimal solution is an optimally scaled lda for maximum  separability between classes that can results
#' in more accurate classification. These direction are refered to as FENN directions.
#' 
#' Specifying the \code{prior} will affect the classification unless over-ridden in \code{predict.fenn}.   
#' 
#' @note This function may be called giving either a formula and optional data frame, or a matrix and
#' grouping factor as the first two arguments.  All other arguments are optional, but \code{subset=} and
#' \code{na.action=}, if required, must be fully named.
#'
#' If a formula is given as the principal argument the object may be modified using \code{update()} in 
#' the usual way.
#' 
#' @seealso \code{\link{predict.fenn}}
#' 
#' @examples
#' Iris <- data.frame(rbind(iris3[,,1], iris3[,,2], iris3[,,3]),
#' Sp = rep(c("s","c","v"), rep(50,3)))
#' train <- sample(1:150, 75)
#' table(Iris$Sp[train])
#' ## your answer may differ
#' ##  c  s  v
#' ## 22 23 30
#' z <- fenn(Sp ~ ., Iris, prior = c(1,1,1)/3, subset = train)
#' predict(z, Iris[-train, ])$class
#' ##  [1] s s s s s s s s s s s s s s s s s s s s s s s s s s s c c c
#' ## [31] c c c c c c c v c c c c v c c c c c c c c c c c c v v v v v
#' ## [61] v v v v v v v v v v v v v v v
#' (z1 <- update(z, . ~ . - Petal.W.))
#' 
#' @export
#' 

fenn <- function(x, ...) UseMethod("fenn")

#' @export
#'
fenn.formula <- function(formula, data, ..., subset, na.action)
{
  m <- match.call(expand.dots = FALSE)
  m$... <- NULL
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  grouping <- model.response(m)
  x <- model.matrix(Terms, m)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)
  if(xint > 0L) x <- x[, -xint, drop = FALSE]
  res <- fenn.default(x, grouping, ...)
  res$terms <- Terms
  ## fix up call to refer to the generic, but leave arg name as `formula'
  cl <- match.call()
  cl[[1L]] <- as.name("fenn")
  res$call <- cl
  res$contrasts <- attr(x, "contrasts")
  res$xlevels <- .getXlevels(Terms, m)
  res$na.action <- attr(m, "na.action")
  res
}

#' @export
#'
fenn.data.frame <- function(x, ...)
{
  res <- fenn(structure(data.matrix(x), class = "matrix"), ...)
  cl <- match.call()
  cl[[1L]] <- as.name("fenn")
  res$call <- cl
  res
}

#' @export
#'
fenn.matrix <- function(x, grouping, ..., subset, na.action)
{
  if(!missing(subset)) {
    x <- x[subset, , drop = FALSE]
    grouping <- grouping[subset]
  }
  if(!missing(na.action)) {
    dfr <- na.action(structure(list(g = grouping, x = x),
                               class = "data.frame"))
    grouping <- dfr$g
    x <- dfr$x
  }
  #    res <- NextMethod("fenn")
  res <- fenn.default(x, grouping, ...)
  cl <- match.call()
  cl[[1L]] <- as.name("fenn")
  res$call <- cl
  res
}

#' @export
#'
model.frame.fenn <- function(formula, ...)
{
  oc <- formula$call
  oc$prior <- oc$tol <- oc$method <- oc$CV <- oc$nu <- NULL
  oc[[1L]] <- quote(stats::model.frame)
  if(length(dots <- list(...))) {
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0L)]
    oc[names(nargs)] <- nargs
  }
  if (is.null(env <- environment(formula$terms))) env <- parent.frame()
  eval(oc, env)
}

#' This is the main function for fitting \code{fenn}.
#'
#' @usage fenn.default(x, ...)
#'                  
#' @param x An object of class \code{"fenn"}.
#' 
#' @details This function is a method for the generic function
#' \code{print()} for class \code{"fenn"}.
#' It can be invoked by calling \code{print(x)} for an
#' object \code{x} of the appropriate class, or directly by
#' calling \code{print.fenn(x)} regardless of the
#' class of the object.
#' @export
#'
fenn.default <- function(x, grouping, prior, tol = 1.0e-4, ...){
  if(is.null(dim(x))) stop("'x' is not a matrix")
  x <- as.matrix(x)
  if(any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  n <- nrow(x)
  p <- ncol(x)
  if(n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- lev1 <- levels(g)
  counts <- as.vector(table(g))
  if(!missing(prior)) {
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != nlevels(g)) stop("'prior' is of incorrect length")
    prior <- prior[counts > 0L]
  }
  if(any(counts == 0L)) {
    empty <- lev[counts == 0L]
    warning(sprintf(ngettext(length(empty),
                             "group %s is empty",
                             "groups %s are empty"),
                    paste(empty, collapse = " ")), domain = NA)
    lev1 <- lev[counts > 0L]
    g <- factor(g, levels = lev1)
    counts <- as.vector(table(g))
  }
  proportions <- counts/n
  if(missing(prior)) {
    prior <- proportions
  }
  
  ng <- length(proportions)
  names(prior) <- names(counts) <- lev1
  
  ## drop attributes to avoid e.g. matrix() methods
  group.means <- tapply(c(x), list(rep(g, p), col(x)), mean)
  means <- colSums(group.means)
  x <- scale(x, center = means, scale = FALSE)
  
  # f1 <- sqrt(diag(var(x - group.means[g,  ])))
  # if(any(f1 < tol)) {
  #   const <- format((1L:p)[f1 < tol])
  #   stop(sprintf(ngettext(length(const),
  #                         "variable %s appears to be constant within groups",
  #                         "variables %s appear to be constant within groups"),
  #                paste(const, collapse = " ")),
  #        domain = NA)
  # }
  
  ## compute scatter matricies
  X.S <- matrix(0, nrow = p, ncol = p)
  X.D <- matrix(0, nrow = p, ncol = p)
  X.C <- matrix(0, nrow = ng, ncol = p)
  for(i in 1:ng){
    X.g <- x[g == lev[i], ]
    X.C[i,] <- colMeans(X.g)
    X.S <- X.S + ((nrow(X.g) - 1) / nrow(X.g)) * cov(X.g)
  }
  
  X.S <- X.S / ng
  
  E <- eigen(X.S,symmetric = T)
  lambda <- E$values
  lambda[which(lambda < tol^2)] <- 0
  
  X.S <- E$vectors %*% diag(lambda) %*% t(E$vectors)
  
  if(ng > 1){
    X.D <- X.S + cov(X.C)
    E <- eigen(X.D,symmetric = T)
    lambda <- E$values
    lambda[which(abs(lambda) < tol^2)] <- 0
    X.D <- E$vectors %*% diag(lambda) %*% t(E$vectors)
  }else{
    X.D <- diag(p)
  }
  
  diag.values <- lambda
  diag.values[which(lambda > 0)] <- 1/sqrt(lambda[which(lambda > 0)])
  X_D_neg_1_2 <- E$vectors %*% diag(diag.values) %*% t(E$vectors)
  
  ## Compute the hamiltonian
  X_tilde_S <- X_D_neg_1_2 %*% X.S %*% X_D_neg_1_2
  E <- eigen(X_tilde_S,symmetric = T)
  lambda <- E$values
  if(ng > 1){
    lambda[which(abs(lambda) < tol^2)] <- 1
  }
  
  X_tilde_S <- E$vectors %*% diag(lambda) %*% t(E$vectors)
  E_tilde_S <- eigen(X_tilde_S,symmetric = T)
  U <- E_tilde_S$vectors
  dd <- E_tilde_S$values
  if(ng > 1){
    dd[which(abs(dd - 1) < tol^2)] <- as.integer(1)
  }
  
  
  ## Temprature vector
  muVec <- -80:80
  #muVec <- -80:(floor(log10(max(lambda))) + 1)
  muVec <- 10^(muVec/10)
  muVec <- sort(unique(c(muVec, dd)))
  
  ## Compute Boltzman distribution
  w <- matrix(0, nrow = length(muVec), ncol = length(dd)) ## each row is a set of weights for the particular mu
  ## Identify infinite values
  l.ind <- which(unlist(lapply(muVec, function(x) any(abs((1/x) * dd) > 100))))
  if(length(l.ind) > 0){
    for(ll in l.ind){
      w[ll, ] <- sapply(dd, function(x) 1 / sum(exp((-1/muVec[ll]) * (dd - x))))
    } 
    if(length(l.ind) < length(muVec)){
      w[-l.ind, ] <- t(sapply(muVec[-l.ind], function(x) exp(((-1/x) * dd)) / sum(exp(((-1/x) * dd)))))
    }
  }else{
    w <- t(as.matrix(sapply(muVec, function(x) exp(((-1/x) * dd)) / sum(exp(((-1/x) * dd))))))
  }
  
  W <- list(w = w, d = dd)  
  
  ## Compute and find the maxima of Fisher Information
  E  <- matrix(W$w, ncol = length(W$d)) %*% matrix(W$d, ncol = 1)
  E2 <- matrix(W$w, ncol = length(W$d)) %*% (matrix(W$d, ncol = 1)^2)
  dE <- (1/(muVec^2))*(E2 - E^2)
  d_ind <- sort(match(unique(W$d), muVec))
  best.mu <- numeric()
  test.mu <- which.max(dE[1:d_ind[1]])
  if (test.mu != d_ind[1]){
    best.mu <- c(test.mu)
  }
  for(i in 1:(length(d_ind)-1)){
    test.mu <- which.max(dE[d_ind[i]:d_ind[i+1]]) + d_ind[i] - 1
    if ((test.mu != d_ind[i]) & (test.mu != d_ind[i+1])){
      best.mu <- c(best.mu, test.mu)
    }
    test.mu <- which.max(dE[(d_ind[i]-1):(d_ind[i]+1)]) + d_ind[i] - 2
    if (test.mu == d_ind[i]){
      best.mu <- c(best.mu, test.mu)
    }
  }
  test.mu <- which.max(dE[(d_ind[length(d_ind)]-1):(d_ind[length(d_ind)]+1)]) + d_ind[length(d_ind)] - 2
  if (test.mu == d_ind[length(d_ind)]){
    best.mu <- c(best.mu, test.mu)
  }
  
  best.mu <- best.mu[length(best.mu)]
  
  ## Compute the final transformation
  diag.vals <- W$w[best.mu,]
  S <- U %*% diag(diag.vals) %*% t(U)
  S1_2 <- U %*% diag(sqrt(diag.vals)) %*% t(U)
  
  scalings <- eigen(S1_2, symmetric = T)
  
  ## Informative Features
  if(ng > 1){
    num.inf.features <- sum(abs(1 - dd) > tol^2)
    informative.dims <- 1:num.inf.features
  }else{
    informative.dims <- 1:which(cumsum(diag.vals[order(diag.vals, decreasing = T)]) > 0.9)[1]
  }
  
  
  ## transform the data
 
  x.tilde <- t(apply(x, 1, function(xx) matrix(X_D_neg_1_2 %*% matrix(xx, nco = 1), nrow = 1))) 
  x.fenn <- x.tilde %*% t(S1_2)
  x.fenn <- t(t(scalings$vectors) %*% t(x.fenn))
  
  
  if(!is.null(dimnames(x))){
    dimnames(group.means)[[2L]] <- colnames(x)
  }
  
  
  
  cl <- match.call()
  cl[[1L]] <- as.name("fenn")
  s <- structure(list(prior = prior, counts = counts, means = group.means,
                      X.S = X.S, X.D = X.D, S1_2 = S1_2, X_D_neg_1_2 = X_D_neg_1_2, informative.dims = informative.dims, 
                      muVec = muVec, best.mu = best.mu, lambdas = W$w, E = E, dE = dE, X_tilde_S = X_tilde_S,
                      x.tilde = x.tilde, x.fenn = x.fenn, grouping = grouping, scalings = scalings, lev = lev,
                      N = n, call = cl),
                 class = "fenn")
  return(s)
}


#' Classify Multivariate Observations 
#'
#' @description Classify multivariate observations in conjunction with \code{fenn}, and also
#' project data onto the scaled linear discriminants.
#'
#' @usage predict.fenn(object, newdata, prior = object$prior, dimen,
#'                               method = c("knn", "lda"), \dots)
#' 
#' @param object Object  of class \code{"fenn"}.
#' @param newdata Data frame of cases to be classified or, if \code{object}
#'                has a formula, a data frame with columns of the same names as the
#'                variables used.  A vector will be interpreted
#'                as a row vector.  If newdata is missing, an attempt will be
#'                made to retrieve the data used to fit the \code{fenn} object.
#' @param reduce.dm The dimension of the data will be reduced before performing knn.
#'
#' @param method This determines how the parameter estimation is handled. With \code{"plug-in"}
#'              (the default) the usual unbiased parameter estimates are used and
#'              assumed to be correct. With \code{"debiased"} an unbiased estimator of
#'              the log posterior probabilities is used, and with \code{"predictive"} the
#'              parameter estimates are integrated out using a vague prior.
#'
#' @param \dots Arguments based from or to other methods
#'
#' @return A list with components:
#' \item{class}{The MAP classification (a factor).}
#' \item{posterior}{If classification method is lda, the posterior probabilities for the 
#'                  classes in the tilde transformed space will be returned.}
#' 
#' @details This function is a method for the generic function \code{predict()} for
#' class \code{"fenn"}.  It can be invoked by calling \code{predict(x)} for
#' an object \code{x} of the appropriate class, or directly by calling
#' \code{predict.fenn(x)} regardless of the class of the object.
#'
#' Missing values in \code{newdata} are handled by returning \code{NA} if the
#' scaled linear discriminants cannot be evaluated. If \code{newdata} is omitted and
#' the \code{na.action} of the fit omitted cases, these will be omitted on the
#' prediction.
#'
#' This version centres the scaled linear discriminants so that the
#' weighted mean (weighted by \code{prior}) of the group centroids is at
#' the origin. 
#' 
#' @examples tr <- sample(1:50, 25)
#' train <- rbind(iris3[tr,,1], iris3[tr,,2], iris3[tr,,3])
#' test <- rbind(iris3[-tr,,1], iris3[-tr,,2], iris3[-tr,,3])
#' cl <- factor(c(rep("s",25), rep("c",25), rep("v",25)))
#' z <- fenn(train, cl)
#' predict(z, test)$class
#' 
#' @export
predict.fenn <- function(object, newdata, method = c('knn', 'lda'), reduce.dim = F, prior = object$prior, ...){
  if(!inherits(object, "fenn")) stop("object not of class \"fenn\"")
  method <- match.arg(method)
  if(!is.null(Terms <- object$terms)) { # formula fit
    Terms <- delete.response(Terms)
    if(missing(newdata)) newdata <- model.frame(object)
    else {
      newdata <- model.frame(Terms, newdata, na.action=na.pass,
                             xlev = object$xlevels)
      if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, newdata)
    }
    x <- model.matrix(Terms, newdata, contrasts = object$contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    if(xint > 0L) x <- x[, -xint, drop = FALSE]
  } else { # matrix or data-frame fit
    if(missing(newdata)) {
      if(!is.null(sub <- object$call$subset))
        newdata <-
          eval.parent(parse(text = paste(deparse(object$call$x,
                                                 backtick = TRUE),
                                         "[", deparse(sub, backtick = TRUE),",]")))
      else newdata <- eval.parent(object$call$x)
      if(!is.null(nas <- object$call$na.action))
        newdata <- eval(call(nas, newdata))
    }
    if(is.null(dim(newdata)))
      dim(newdata) <- c(1L, length(newdata))  # a row vector
    x <- as.matrix(newdata)		# to cope with dataframes
  }
  
  if(ncol(x) != ncol(object$means)) stop("wrong number of variables")
  if(length(colnames(x)) > 0L &&
     any(colnames(x) != dimnames(object$means)[[2L]]))
    warning("variable names in 'newdata' do not match those in 'object'")
  ng <- length(object$prior)
  if(!missing(prior)) {
    if(any(prior < 0) || round(sum(prior), 5) != 1) stop("invalid 'prior'")
    if(length(prior) != ng) stop("'prior' is of incorrect length")
  }
  ## remove overall means to keep distances small
  #means <- colSums(prior*object$means)
  means <- colSums(object$means)
  x <- scale(x, center = means, scale = FALSE)
  
  S1_2 <- object$S1_2
  X_D_neg_1_2 <- object$X_D_neg_1_2
  scalings <- object$scalings
  
  x.tilde <- t(apply(x, 1, function(xx) matrix(X_D_neg_1_2 %*% matrix(xx, nco = 1), nrow = 1))) 
  x.fenn <- x.tilde %*% t(S1_2)
  x.fenn <- t(t(scalings$vectors) %*% t(x.fenn))
  
  
  if(method == "knn"){
    if(reduce.dim){
      knn.pred <- knn(object$x.fenn[,object$informative.dims], x.fenn[,object$informative.dims], 
                      factor(object$grouping), k = 3, prob=F)
    }else{
      knn.pred <- knn(object$x.fenn, x.fenn, factor(object$grouping), k = 3, prob=F)
    }
    class.label <- knn.pred
    posterior <- NULL
  }else if(method == 'lda'){
    D.train <- data.frame(y = object$grouping, object$x.tilde)
    fit.lda <- lda(y~., data = D.train)
    D.test <- data.frame(x.tilde)
    predict.lda <- predict(fit.lda, newdata = D.test)
    class.label <- predict.lda$class
    posterior <- predict.lda$posterior
  }
  #dimnames(posterior) <- rownames(x)
  
  list(class = class.label, posterior = posterior, x.fenn = x.fenn)
}

#' Print Method for Class \code{fenn}.
#' 
#' Prints summary table of fenn fitted object.
#'
#' @usage print.fenn(x, ...)
#'                  
#' @param x An object of class \code{"fenn"}.
#' 
#' @details This function is a method for the generic function
#' \code{print()} for class \code{"fenn"}.
#' It can be invoked by calling \code{print(x)} for an
#' object \code{x} of the appropriate class, or directly by
#' calling \code{print.fenn(x)} regardless of the
#' class of the object.
#' @export
#'
print.fenn <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    names(cl)[2L] <- ""
    cat("Call:\n")
    dput(cl, control = NULL)
  }
  cat("\nPrior probabilities of groups:\n")
  print(x$prior, ...)
  cat("\nGroup means:\n")
  print(x$means, ...)
  cat("\nCoefficients of linear discriminants:\n")
  print(x$scalings$values, ...)
  cat("\nInformative dimentions:\n")
  print(x$informative.dims, ...)
  invisible(x)
}


#' Plot Method for Class \code{fenn}.
#' 
#' Plots a set of data on one, two or more scaled linear discriminants.
#'
#' @usage plot.fenn(x, panel = panel.fenn, \dots, cex = 0.7, dimen,
#'                  abbrev = FALSE, xlab = "FENN1", ylab = "FENN2")
#'                  
#' @param x An object of class \code{"fenn"}.
#' 
#' @param panel The panel function used to plot the data.
#'
#' @param \dots Additional arguments to \code{pairs}, \code{ldahist} or \code{eqscplot}.
#' 
#' @param cex Graphics parameter \code{cex} for labels on plots.
#'
#' @param abbrev Whether the group labels are abbreviated on the plots. If \code{abbrev > 0}
#'               this gives \code{minlength} in the call to \code{abbreviate}.
#'               
#' @param xlab Label for the \code{x} axis.
#'
#' @param ylab Label for the \code{y} axis.
#' 
#' @details This function is a method for the generic function
#' \code{plot()} for class \code{"fenn"}.
#' It can be invoked by calling \code{plot(x)} for an
#' object \code{x} of the appropriate class, or directly by
#' calling \code{plot.fenn(x)} regardless of the
#' class of the object.
#'
#' The behaviour is determined by the value of \code{dimen}. For
#' \code{dimen > 2}, a \code{pairs} plot is used. For \code{dimen = 2}, an
#' equiscaled scatter plot is drawn. For \code{dimen = 1}, a set of
#' histograms or density plots are drawn.  Use argument \code{type} to
#' match \code{"histogram"} or \code{"density"} or \code{"both"}.
#'
#' @export
#'
plot.fenn <- function(x, panel = panel.fenn, ..., cex = 0.7,
                      dimen, abbrev = FALSE,
                      xlab = "FENN1", ylab = "FENN2", reduce.dm = T)
{
  panel.fenn <- function(x, y, ...) text(x, y, as.character(g), cex = cex, ...)
  if(!is.null(Terms <- x$terms)) { # formula fit
    data <- model.frame(x)
    g <- model.response(data)
  } else { # matrix or data-frame fit
    xname <- x$call$x
    gname <- x$call[[3L]]
    if(!is.null(sub <- x$call$subset)) {
      g <- eval.parent(parse(text=paste(deparse(gname, backtick=TRUE),
                                        "[", deparse(sub, backtick=TRUE),"]")))
    } else {
      g <- eval.parent(gname)
    }
  }
  if(abbrev) levels(g) <- abbreviate(levels(g), abbrev)
  means <- colMeans(x$means)
  X <- x$x.fenn
  colnames(X) <- paste('FENN', 1:ncol(X), sep = '')
  if(reduce.dm){
    X <- X[,x$informative.dims, drop = F]
  }
  colors = rainbow(length(unique(g)))
  names(colors) = unique(g)
  
  if(ncol(X) > 2L) {
    #pairs(X, panel = panel, col = colors[g], ...)
    pairs(X, col = colors[g],pch = 20, ...)
  } else if(ncol(X) == 2L)  {
    eqscplot(X[, 1L:2L], xlab = xlab, ylab = ylab, type = "n",  ...)
    panel(X[, 1L], X[, 2L], col = colors[g],pch = 20, ...)
  } 
  invisible(NULL)
}
