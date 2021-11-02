
ccr_smreg <- function(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, var = c("None", "CF", "BS"), nboot = 1, ...){
  UseMethod("ccr_smreg")
}

ccr_smreg.default <- function(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, var = c("None", "BS", "CF"), nboot = 1, ...){
  nc <- length(unique(data$clusterid))
  x <-  sort(unique(data$x))
  t <- sort(unique(data$x[data$c>0 & data$r==1]))
  Q1 <- as.numeric(quantile(data$x[data$c > 0], 0.05))
  Q9 <- as.numeric(quantile(data$x[data$c > 0], 0.95))
  q1 <- as.numeric(quantile(t, 0.05))
  q9 <- as.numeric(quantile(t, 0.95))

  est1 <- clustered_mpple_est(data, formula11, formula21, cause = 1, w = w, t = t)
  est2 <- clustered_mpple_est(data, formula12, formula22, cause = 2, w = w, t = t)

  var <- match.arg(var, c("None", "CF", "BS"))

  if(var == "None"){
    se1 <- NA
    se2 <- NA
  } else if (var == "BS"){
    message("please wait...")
    se <- clustered_mpple_se_bs(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = w, nboot = nboot, x = x, t = t)
    se1 <- se[[1]]
    se2 <- se[[2]]
  } else if (var == "CF"){
    message("please wait...")
    se1 <- clustered_mpple_se_cf(data, formula11, formula21, cause = 1, w = w)
    se2 <- clustered_mpple_se_cf(data, formula12, formula22, cause = 2, w = w)
  }
  formula <- list(formula11, formula12, formula21, formula22)
  out <- list(est1 = est1, est2 = est2, se1 = se1, se2 = se2, formula = formula, x = x, t = t, Q1 = Q1, Q9 = Q9, q1 = q1, q9 = q9, nc = nc)
  class(out) <- "ccr_smreg"
  attr(out, "var") <- var
  return(out)
}

summary.ccr_smreg <- function(x,...){

  var <- attr(x, "var")
  var <- match.arg(var, c("None", "CF", "BS"))

  beta1 <- x$est1$beta
  hr1 <- exp(beta1)
  beta2 <- x$est2$beta
  hr2 <- exp(beta2)

  if (var == "None"){
    cause1 <- cbind(beta1, hr1)
    row.names(cause1) <- all.vars(formula[[3]])[-c(1,2)]
    colnames(cause1) <- c("beta", "exp(beta)")
    cause2 <- cbind(beta2, hr2)
    row.names(cause2) <- all.vars(formula[[4]])[-c(1,2)]
    colnames(cause2) <- c("beta", "exp(beta)")
    out <- list(cause1, cause2)
  } else if (var == "CF"){
    se1 <- x$se1$SE_B
    se2 <- x$se2$SE_B
    CI_lower1 <- exp(beta1 - 1.96*se1)
    CI_lower2 <- exp(beta2 - 1.96*se2)
    CI_upper1 <- exp(beta1 + 1.96*se1)
    CI_upper2 <- exp(beta2 + 1.96*se2)
    cause1 <- cbind(beta1, se1, hr1, CI_lower1, CI_upper1)
    row.names(cause1) <- all.vars(formula[[3]])[-c(1,2)]
    colnames(cause1) <- c("beta", "se","exp(beta)", "95% CI lower limit", "95% CI upper limit")
    cause2 <- cbind(beta2, se2, hr2, CI_lower2, CI_upper2)
    row.names(cause2) <- all.vars(formula[[4]])[-c(1,2)]
    colnames(cause2) <- c("beta", "se","exp(beta)", "95% CI lower limit", "95% CI upper limit")
    out <- list(cause1, cause2)
  } else if (var == "BS"){
    se1 <- apply(x$se1$beta_boot, 2, sd)
    se2 <- apply(x$se2$beta_boot, 2, sd)
    CI_lower1 <- exp(beta1 - 1.96*se1)
    CI_lower2 <- exp(beta2 - 1.96*se2)
    CI_upper1 <- exp(beta1 + 1.96*se1)
    CI_upper2 <- exp(beta2 + 1.96*se2)
    cause1 <- cbind(beta1, se1, hr1, CI_lower1, CI_upper1)
    row.names(cause1) <- all.vars(formula[[3]])[-c(1,2)]
    colnames(cause1) <- c("beta", "se","exp(beta)", "95% CI lower limit", "95% CI upper limit")
    cause2 <- cbind(beta2, se2, hr2, CI_lower2, CI_upper2)
    row.names(cause2) <- all.vars(formula[[4]])[-c(1,2)]
    colnames(cause2) <- c("beta", "se","exp(beta)", "95% CI lower limit", "95% CI upper limit")
    out <- list(cause1, cause2)
  }

  class(out) <- "summary.ccr_smreg"
  return(out)
}

print.summary.ccr_smreg <- function(x,...){
  cat("\n")
  cat("cause 1: \n")
  print(out[[1]])
  cat("\n")
  cat("cause 2: \n")
  print(out[[2]])
}
