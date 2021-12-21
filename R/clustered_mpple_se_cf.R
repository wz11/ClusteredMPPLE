


#' Obtain coefficient estimates standard error using close form variance
#'
#' @description Obtain coefficient estimates standard error for semiparametric regression for competing risks data
#' with missing cause of failure using close form variance.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param data a data frame in suitable format.
#' @param formula1 a formula relating the response for \code{y} to a set of covariates.
#' @param formula2 a formula relating the survival object \code{Surv(x, d)} to a set of covariates.
#' @param cause 1 or 2, the cause of failure to be evaluated.
#' @param w logical value: if TRUE, the estimation procedure is weighed by cluster size.
#' @param ... for future methods.
#'
#' @return A list with components:
#' \item{beta}{a vector of regression coefficients}
#' \item{SE_B}{a vector of standard errors of regression coefficients}
#' \item{time}{a vector of time}
#' \item{H}{a vector of cumulative hazard function}
#' \item{phi}{useful quantities for incluence function}
#' \item{dphi}{useful quantities for incluence function}
#'
clustered_mpple_se_cf <- function(data, formula1 =  y ~ x + Z1 + Z2, formula2 = Surv(x, d) ~ Z1 + Z2, cause = 1, w = TRUE, ...){

  n <- dim(data)[1]
  nc <- length(unique(data$clusterid))

  data$include <- 1*(data$r==1 & data$c>0)/(1*(w==FALSE)+data$clustersize*(w==TRUE))
  data$y <- 1*(data$c==cause)

  model <- geepack::geeglm(formula1, family = "binomial", data = data, id = clusterid, weights = include, corstr = "independence")
  data$yhat <- stats::predict(model, data, type = "response")
  omega <- t(model$geese$infls)[ ,-(dim(model$geese$infls)[1])] * nc

  # weighted cox proportional hazard model
  data$d <- data$r*(data$c == cause) + (1-data$r)*(data$yhat > 0)
  data$weight <- data$r + (1-data$r)*data$yhat
  data$weight <- data$weight + (data$weight == 0)

  dt0 <- data[data$r==0, ]
  dt0$weight <- 1 - dt0$weight
  dt0$d <- 0

  data1 <- rbind(data, dt0)
  data1$weight <- data1$weight/(1*(w==FALSE)+data1$clustersize*(w==TRUE))

  mod <- survival::coxph(formula = formula2, weights = weight, data = data1)
  beta <- stats::coef(mod)
  Haz <- survival::basehaz(mod, centered = FALSE)
  time <- Haz$time
  H <- Haz$hazard

  H_b <- stats::vcov(mod)*n
  data$w <- data$d*data$weight
  cov <- all.vars(formula2)[-c(1,2)]

  ar_t <- function(t){
    sum((data$x >= t) * exp(as.matrix(data[ ,cov]) %*% beta)/(1*(w==FALSE)+data$clustersize*(w==TRUE)))
  }

  E_t <- function(t){
    colSums(data[ ,cov]*(data$x>=t)*exp(as.matrix(data[ ,cov]) %*% beta)/(1*(w==FALSE)+data$clustersize*(w==TRUE)))/ar_t(t)
  }

  dM_t <- function(t){
    data$w*(data$x==t)-(data$x>=t)*exp(as.matrix(data[ ,cov]) %*% beta)*sum(data$w * (data$x==t)/(1*(w==FALSE)+data$clustersize*(w==TRUE)))/ar_t(t)
  }

  # Actual calculations

  E <- sapply(data$x, E_t, simplify = TRUE)
  dM <- sapply(data$x, dM_t, simplify = TRUE) #Each column corresponds to a t

  psi_n <- sapply(1:n, FUN=function(x){colSums(t((as.vector(t(data[x,cov]))-E))*dM[x,])})

  psi_j <- sapply(1:length(cov),FUN=function(x){as.vector(tapply(psi_n[x,], data$clusterid, mean))})


  data$one <- 1
  logcovs <- c("one",  all.vars(formula1)[-1])
  COV <- data[ ,logcovs]

  if (length(cov)==1){
    R_ij <-  as.matrix((1 - data$r)*(data[ ,cov] - t(E))*(data$c > 0)) %*% as.matrix(data$yhat*(1 - data$yhat)*COV/(1*(w==FALSE)+data$clustersize*(w==TRUE)))
  }else{
    R_ij <-  t((1 - data$r)*(data[ ,cov] - t(E))*(data$c > 0)) %*% as.matrix(data$yhat*(1 - data$yhat)*COV/(1*(w==FALSE)+data$clustersize*(w==TRUE)))
  }

  R_j <- t(R_ij/nc)
  psi <- (psi_j + omega %*% R_j) %*% H_b

  m_psi2 <- t(psi) %*% psi/ nc
  V_B <- sqrt(diag(m_psi2 / nc))

  dH_t <- function(t){
    sum(data$w * (data$x == t)/(1*(w==FALSE)+data$clustersize*(w==TRUE))) / ar_t(t)
  }

  dH <- sapply(data$x, dH_t, simplify = TRUE)
  ar <- sapply(data$x, ar_t, simplify = TRUE)/nc

  phi_t <- function(t){
    R_ij <- (1 - data$r)*((data$x <= t)/ar) * (data$c > 0) * (data$yhat * (1 - data$yhat) * COV)/(1*(w==FALSE)+data$clustersize*(w==TRUE))
    R_j <- colSums(R_ij)/nc
    out <- as.vector(tapply(colSums(t(dM[ ,(data$x <= t)]) / ar[data$x <= t]), data$clusterid, mean)) - as.vector(psi %*% as.matrix(colSums(matrix((data$x <= t) * t(E) * dH, ncol = length(cov), nrow = n)))) +
      as.vector(omega %*% R_j)
    return(out)
  }

  phi <- t(sapply(time, phi_t, simplify = TRUE))

  dphi_t <- function(t){
    R_ij <- (1 - data$r)*((data$x == t)/(ar_t(t)/nc)) * (data$c > 0) * (data$yhat * (1 - data$yhat) * COV)/(1*(w==FALSE)+data$clustersize*(w==TRUE))
    R_j <- colSums(R_ij)/nc
    out <- as.vector(tapply(dM_t(t)/(ar_t(t)/nc), data$clusterid, mean)) - psi %*% (as.matrix(E_t(t)) * dH_t(t)) +
      omega %*% R_j
    return(out)
  }
  dphi <- t(sapply(time, dphi_t, simplify = TRUE))

  out <- list(beta = beta, SE_B = V_B, time = time, H = H, phi = phi, dphi = dphi)
  return(out)
}
