

#' Obtain coefficient estimates standard error using bootstrap method
#'
#' @description Obtain coefficient estimates standard error for semiparametric regression for competing risks data
#' with missing cause of failure using bootstrap.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param data a data frame in suitable format.
#' @param formula1 a formula relating the response for cause 1 \code{y} to a set of covariates.
#' @param formula2 a formula relating the survival object \code{Surv(x, d)} to a set of covariates.
#' @param w logical value: if TRUE, the estimation procedure is weighed by cluster size.
#' @param nboot a number of bootstrap samples for estimating variances when bootstrapping methods are used.
#' @param x a vector of time points for cumulative incidence function.
#' @param t a vector of time points for cumulative residual process.
#' @param ... for future methods.
#'
#' @return A list with components:
#' \item{se1}{a list of bootstrap samples for cause 1}
#' \item{se2}{a list of bootstrap samples for cause 1}
#'

clustered_mpple_se_bs <- function(data, formula1 =  y ~ x + Z1 + Z2, formula2 = Surv(x, d) ~ Z1 + Z2, w = TRUE, nboot = 1000, x, t,...){

  b1 <- length(all.vars(formula2)) - 2
  b2 <- length(all.vars(formula2)) - 2
  clusterid <- as.vector(unique(data$clusterid))
  nc <- length(unique(data$clusterid))

  beta_boot_cause1 <- matrix(NA, nboot, b1)
  CIF_boot_cause1 <- matrix(NA, nboot, length(x))
  crp_boot_cause1 <- matrix(NA, nboot, length(t))

  beta_boot_cause2 <- matrix(NA, nboot, b2)
  CIF_boot_cause2 <- matrix(NA, nboot, length(x))
  crp_boot_cause2 <- matrix(NA, nboot, length(t))

  for (i in 1:nboot){
    print(i)
    set.seed(i)
    cluster <- sample(clusterid, nc, replace = TRUE)
    data_temp <- data.table::rbindlist(lapply(cluster, function(x) data[data$clusterid == x,]))

    result_temp1 <- clustered_mpple_est(data = data_temp, formula1 =  formula1, formula2 = formula2, cause = 1, w = w, t = t)
    result_temp2 <- clustered_mpple_est(data = data_temp, formula1 =  formula1, formula2 = formula2, cause = 2, w = w, t = t)

    cif_temp1 <- clustered_cif_est(est1 = result_temp1, est2 = result_temp2, cause = 1, t=x)$CIF
    cif_temp2 <- clustered_cif_est(est1 = result_temp1, est2 = result_temp2, cause = 2, t=x)$CIF

    beta_boot_cause1[i,] <- result_temp1$beta
    beta_boot_cause2[i,] <- result_temp2$beta

    crp_boot_cause1[i,] <- result_temp1$crp
    crp_boot_cause2[i,] <- result_temp2$crp

    CIF_boot_cause1[i,] <- cif_temp1
    CIF_boot_cause2[i,] <- cif_temp2

  }

  se1 <- list(beta_boot = beta_boot_cause1, crp_boot = crp_boot_cause1, CIF_boot = CIF_boot_cause1)
  se2 <- list(beta_boot = beta_boot_cause2, crp_boot = crp_boot_cause2, CIF_boot = CIF_boot_cause2)

  out <- list(se1, se2)
  return(out)
}
