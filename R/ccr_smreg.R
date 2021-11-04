
#' Semiparametric Marginal Regression for Clustered Competing Risks Data with Missing Cause of Failure
#'
#' @description This function \code{ccr_smreg} conducts semiparametric marginal regression for clustered competing risks data with missing cause of failure. It implements the proposed maximum partial pseudolikelihood estimation method by Zhou, W. et al.(2021) based on semiparametric marginal proportional cause-specific hazards model. The standard errors for the estimated regression coefficients are calculated based on the close-form variance or bootstrapping.

#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param data a data frame in suitable format.
#' @param formula11 a formula relating the response for cause 1 \code{y} to a set of covariates.
#' @param formula12 a formula relating the response for cause 1 \code{y} to a set of covariates.
#' @param formula21 a formula relating the survival object \code{Surv(x, d)} to a set of covariates.
#' @param formula22 a formula relating the survival object \code{Surv(x, d)} to a set of covariates.
#' @param w logical value: if TRUE, the estimation procedure is weighed by cluster size.
#' @param var.method a character string specifying the method for variance calcultaion. If None, the variance will
#' not be calculated; If CF, the close form variance will be used; If BS, the bootstrapping methods will be
#' used.
#' @param nboot a number of bootstrap samples for estimating variances when bootstrapping methods are used
#' @param ... for future methods
#'
#' @return An object of class \code{ccr_smreg} with components:
#' \item{est1}{a list for the estimates for cause 1}
#' \item{est2}{a list for the estimates for cause 2}
#' \item{se1}{a list for the standard errors for cause 1}
#' \item{se2}{a list for the standard errors for cause 2}
#' \item{formula}{a list for the formulas}
#' \item{x}{a vector for time points where cumulative incidence functions are evaluated}
#' \item{t}{a vector for time points where cumulative residual processes are evaluated}
#' \item{Q1}{lower limit for cumulative incidence functions}
#' \item{Q9}{upper limit for cumulative incidence functions}
#' \item{q1}{lower limit for cumulative residual process}
#' \item{q9}{upper limit for cumulative residual process}
#' \item{nc}{number of clusters}
#'
#' @references
#' {Zhou, W., Bakoyannis, G., Zhang, Y., & Yiannoutsos, C. T. (2021). Semiparametric Marginal Regression for
#' Clustered Competing Risks Data with Missing Cause of Failure. arXiv preprint arXiv:2104.09090.}
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(ClusteredMPPLE)
#' data("simulated_data")
#' fit <- ccr_smreg(data=simulated_data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2,
#' formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2,
#' w = TRUE, var.method = "BS", nboot = 10)
#' summary(fit)
#' }
#' @export
ccr_smreg <- function(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, var.method = c("None", "CF", "BS"), nboot = 1, ...){
  UseMethod("ccr_smreg")
}

#' @export
ccr_smreg.default <- function(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, var.method = c("None", "BS", "CF"), nboot = 1, ...){
  nc <- length(unique(data$clusterid))
  x <-  sort(unique(data$x))
  t <- sort(unique(data$x[data$c>0 & data$r==1]))
  Q1 <- as.numeric(stats::quantile(data$x[data$c > 0], 0.1))
  Q9 <- as.numeric(stats::quantile(data$x[data$c > 0], 0.9))
  q1 <- as.numeric(stats::quantile(t, 0.1))
  q9 <- as.numeric(stats::quantile(t, 0.9))

  est1 <- clustered_mpple_est(data, formula11, formula21, cause = 1, w = w, t = t)
  est2 <- clustered_mpple_est(data, formula12, formula22, cause = 2, w = w, t = t)

  var.method <- match.arg(var.method, c("None", "CF", "BS"))

  if(var.method == "None"){
    se1 <- NA
    se2 <- NA
  } else if (var.method == "BS"){
    message("please wait...")
    se <- clustered_mpple_se_bs(data, formula11, formula12, formula21, formula22, w = w, nboot = nboot, x = x, t = t)
    se1 <- se[[1]]
    se2 <- se[[2]]
  } else if (var.method == "CF"){
    message("please wait...")
    se1 <- clustered_mpple_se_cf(data, formula11, formula21, cause = 1, w = w)
    se2 <- clustered_mpple_se_cf(data, formula12, formula22, cause = 2, w = w)
  }
  formula <- list(formula11 = formula11, formula12 = formula12, formula21 = formula21, formula22 = formula22)
  out <- list(est1 = est1, est2 = est2, se1 = se1, se2 = se2, formula = formula, x = x, t = t, Q1 = Q1, Q9 = Q9, q1 = q1, q9 = q9, nc = nc)
  class(out) <- "ccr_smreg"
  attr(out, "var.method") <- var.method
  return(out)
}





#' Summary method for \code{ccr_smreg} object
#'
#' @description Produce a summary of the fitted semiparametric marginal regression model
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param object the object of a fitted ccr_smreg model
#' @param ... for future methods
#'
#' @return An object of class \code{summary.ccr_smreg} with components:
#' \item{cause1}{a dataframe for the model summary for cause 1}
#' \item{cause2}{a dataframe for the model summary for cause 2}
#'
#' @seealso \code{\link[ClusteredMPPLE]{print.summary.ccr_smreg}} for printing the summarized results
#'
#' @examples
#' \dontrun{
#' summary(fit)
#' }
#' @export
summary.ccr_smreg <- function(object,...){

  var.method <- attr(object, "var.method")
  var.method <- match.arg(var.method, c("None", "CF", "BS"))

  beta1 <- object$est1$beta
  hr1 <- exp(beta1)
  beta2 <- object$est2$beta
  hr2 <- exp(beta2)

  if (var.method == "None"){
    cause1 <- cbind(beta1, hr1)
    row.names(cause1) <- all.vars(object$formula$formula21)[-c(1,2)]
    colnames(cause1) <- c("beta", "exp(beta)")
    cause2 <- cbind(beta2, hr2)
    row.names(cause2) <- all.vars(object$formula$formula22)[-c(1,2)]
    colnames(cause2) <- c("beta", "exp(beta)")
    out <- list(cause1, cause2)
  } else if (var.method == "CF"){
    se1 <- object$se1$SE_B
    se2 <- object$se2$SE_B
    p.value1 <- 2*stats::pnorm(-abs(beta1/se1))
    p.value2 <- 2*stats::pnorm(-abs(beta2/se2))
    CI_lower1 <- exp(beta1 - 1.96*se1)
    CI_lower2 <- exp(beta2 - 1.96*se2)
    CI_upper1 <- exp(beta1 + 1.96*se1)
    CI_upper2 <- exp(beta2 + 1.96*se2)
    cause1 <- cbind(beta1, se1, p.value1, hr1, CI_lower1, CI_upper1)
    row.names(cause1) <- all.vars(object$formula$formula21)[-c(1,2)]
    colnames(cause1) <- c("beta", "se", "p-value", "exp(beta)", "95% CI lower limit", "95% CI upper limit")
    cause2 <- cbind(beta2, se2, p.value2, hr2, CI_lower2, CI_upper2)
    row.names(cause2) <- all.vars(object$formula$formula22)[-c(1,2)]
    colnames(cause2) <- c("beta", "se", "p-value", "exp(beta)", "95% CI lower limit", "95% CI upper limit")
    out <- list(cause1 = cause1, cause2 = cause2)
  } else if (var.method == "BS"){
    se1 <- apply(object$se1$beta_boot, 2, stats::sd)
    se2 <- apply(object$se2$beta_boot, 2, stats::sd)
    p.value1 <- 2*stats::pnorm(-abs(beta1/se1))
    p.value2 <- 2*stats::pnorm(-abs(beta2/se2))
    CI_lower1 <- exp(beta1 - 1.96*se1)
    CI_lower2 <- exp(beta2 - 1.96*se2)
    CI_upper1 <- exp(beta1 + 1.96*se1)
    CI_upper2 <- exp(beta2 + 1.96*se2)
    cause1 <- cbind(beta1, se1, p.value1, hr1, CI_lower1, CI_upper1)
    row.names(cause1) <- all.vars(object$formula$formula21)[-c(1,2)]
    colnames(cause1) <- c("beta", "se", "p-value", "exp(beta)", "95% CI lower limit", "95% CI upper limit")
    cause2 <- cbind(beta2, se2, p.value2, hr2, CI_lower2, CI_upper2)
    row.names(cause2) <- all.vars(object$formula$formula22)[-c(1,2)]
    colnames(cause2) <- c("beta", "se", "p-value", "exp(beta)", "95% CI lower limit", "95% CI upper limit")
    out <- list(cause1 = cause1, cause2 = cause2)
  }

  class(out) <- "summary.ccr_smreg"
  return(out)
}




#' Print method for \code{summary.ccr_smreg} objects
#'
#' @description Produce a printed summary of the fitted semiparametric marginal regression model
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param x the result of a call to 'summary.ccr_smreg'
#' @param ... for future methods
#'
#' @export
print.summary.ccr_smreg <- function(x,...){
  cat("\n")
  cat("cause 1: \n")
  print(x$cause1)
  cat("\n")
  cat("cause 2: \n")
  print(x$cause2)
}
