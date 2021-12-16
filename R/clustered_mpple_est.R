
#' Obtain Coefficient Estimates with proposed method using proposed method
#'
#' @description Obtain coefficient estimates for semiparametric marginal regression for clustered competing risks data with missing cause of failure.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param data a data frame in suitable format.
#' @param formula1 a formula relating the response for \code{y} to a set of covariates.
#' @param formula2 a formula relating the survival object \code{Surv(x, d)} to a set of covariates.
#' @param cause 1 or 2, the cause of failure to be evaluated.
#' @param w logical value: if TRUE, the estimation procedure is weighed by cluster size.
#' @param t a vector of time points.
#' @param ... for future methods.
#'
#' @return A list with components:
#' \item{beta}{a vector of regression coefficients}
#' \item{time}{a vector of time}
#' \item{H}{a vector of cumulative hazard function}
#' \item{crp}{a vector of cumulative residual process}
#'


clustered_mpple_est <- function(data, formula1 =  y ~ x + Z1 + Z2, formula2 = Surv(x, d) ~ Z1 + Z2, cause = 1, w = TRUE, t,...){

  n <- dim(data)[1]
  nc <- length(unique(data$clusterid))


  data$include <- 1*(data$r==1 & data$c>0)/(1*(w==FALSE)+data$clustersize*(w==TRUE))
  data$y <- 1*(data$c==cause)

  model <- stats::glm(formula1, family = "binomial", data = data, weights = include)
  data$yhat <- stats::predict(model, data, type = "response")

  # weighted cox proportional hazard model
  data$d <- data$r*(data$c == cause) + (1-data$r)*(data$yhat > 0)
  data$weight <- data$r + (1-data$r)*data$yhat
  data$weight <- data$weight + (data$weight == 0)

  dt0 <- data[data$r==0, ]
  dt0$weight <- 1 - dt0$weight
  dt0$d <- 0

  data1 <<- rbind(data, dt0)
  data1$weight <- data1$weight/(1*(w==FALSE)+data1$clustersize*(w==TRUE))

  mod <- survival::coxph(formula = formula2, weights = weight, data = data1)
  beta <- stats::coef(mod)
  Haz <- survival::basehaz(mod, centered = FALSE)
  time <- Haz$time
  H <- Haz$hazard

  cum_residual_process <- function(t){
    sum(data$r*((data$x == t)*(data$c == cause) - data$yhat*(data$x == t)*(data$c > 0)) / (1*(w==FALSE)+data$clustersize*(w==TRUE))) / nc
  }

  crp <- sapply(t, cum_residual_process, simplify = TRUE)

  out <- list(beta = beta, time = time, H = H, crp = crp)
  return(out)
}

