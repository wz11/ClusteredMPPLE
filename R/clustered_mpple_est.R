
library(geepack)
library(survival)


clustered_mpple_est <- function(data, formula1 =  y ~ x + Z1 + Z2, formula2 = Surv(x, d) ~ Z1 + Z2, cause = 1, w = TRUE, t, ...){

  n <- dim(data)[1]
  nc <- length(unique(data$clusterid))

  data$include <- 1*(data$r==1 & data$c>0)/(1*(w==FALSE)+data$clustersize*(w==TRUE))
  data$y <- 1*(data$c==cause)

  model <- geepack::geeglm(formula1, family = "binomial", data = data, id = clusterid, weights = include, corstr = "independence")
  data$yhat <- predict(model, data, type = "response")

  # weighted cox proportional hazard model
  data$d <- data$r*(data$c == cause) + (1-data$r)*(data$yhat > 0)
  data$weight <- data$r + (1-data$r)*data$yhat
  data$weight <- data$weight + (data$weight == 0)

  dt0 <- data[data$r==0, ]
  dt0$weight <- 1 - dt0$weight
  dt0$d <- 0

  data1 <<- rbind(data, dt0)
  data1$weight <<- data1$weight/(1*(w==FALSE)+data1$clustersize*(w==TRUE))

  mod <- survival::coxph(formula = formula2, weights = weight, data = data1)
  beta <- coef(mod)
  Haz <- basehaz(mod, centered = FALSE)
  time <- Haz$time
  H <- Haz$hazard

  cum_residual_process <- function(t){
    sum(data$r*((data$x == t)*(data$c == cause) - data$yhat*(data$x == t)*(data$c > 0)) / (1*(w==FALSE)+data$clustersize*(w==TRUE))) / nc
  }

  crp <- sapply(t, cum_residual_process, simplify = TRUE)

  out <- list(beta = beta, time = time, H = H, crp = crp)
  return(out)
}


naive_mpple_est <- function(data, formula1 =  y ~ x + Z1 + Z2, formula2 = Surv(x, d) ~ Z1 + Z2, cause = 1, w = FALSE, t, ...){

  n <- dim(data)[1]
  nc <- length(unique(data$clusterid))

  data$include <- 1*(data$r==1 & data$c>0)/(1*(w==FALSE) + data$clustersize*(w==TRUE))
  data$y <- 1*(data$c==cause)

  model <- glm(formula1, family = "binomial", data = data, weights = include)
  data$yhat <- predict(model, data, type = "response")

  # weighted cox proportional hazard model
  data$d <- data$r*(data$c == cause) + (1-data$r)*(data$yhat > 0)
  data$weight <- data$r + (1-data$r)*data$yhat
  data$weight <- data$weight + (data$weight == 0)

  dt0 <- data[data$r==0, ]
  dt0$weight <- 1 - dt0$weight
  dt0$d <- 0

  data1 <<- rbind(data, dt0)
  data1$weight <<- data1$weight/(1*(w==FALSE) + data1$clustersize*(w==TRUE))

  mod <- survival::coxph(formula = formula2, weights = weight, data = data1)
  beta <- coef(mod)
  Haz <- basehaz(mod, centered = FALSE)
  time <- Haz$time
  H <- Haz$hazard

  cum_residual_process <- function(t){
    sum(data$r*((data$x == t)*(data$c == cause) - data$yhat*(data$x == t)*(data$c > 0)) / (1*(w==FALSE)+data$clustersize*(w==TRUE))) / nc
  }

  crp <- sapply(t, cum_residual_process, simplify = TRUE)

  out <- list(beta = beta, time = time, H = H, crp = crp)
  return(out)
}
