classical.ve <- function(reg.res) { # equivalent to vcov(.)
  X <- model.matrix(reg.res)
  XtX <- crossprod(X)
  A <- solve(XtX)
  B <- XtX * (sum(reg.res$residual^2)/reg.res$df.residual)
  
  A %*% B %*% A
}

hrve <- function(reg.res) { # equivalent to vcovHC(., type="HC0")
  X <- model.matrix(reg.res)
  XtX <- crossprod(X)
  A <- solve(XtX)
  B <- aaply(reg.res$residual * X, 1, function(x) crossprod(t(x))) %>%
    aaply(2:3, sum)
  
  A %*% B %*% A
}

crve <- function(reg.res, clusters) { # equivalent to plm::vcovHC(., type="HC0", cluster="group")
  X <- model.matrix(reg.res)
  XtX <- crossprod(X)
  A <- solve(XtX)
  
  if (missing(clusters)) {
    if (!is(reg.res, "plm")) {
      stop("clusters must be explicity specified or reg.res needs to be a plm object.")
    }
    
    clusters <- index(reg.res)[, 1] 
  }
  
  index.mask <- tapply(seq_len(nrow(X)), clusters)
  
  B <- foreach(index=seq_len(max(index.mask))) %do% {
    X[index.mask == index, ] %>% crossprod(reg.res$residuals[index.mask == index]) %>% t %>% crossprod
  } %>% Reduce("+", .)
  
  A %*% B %*% A
}
  