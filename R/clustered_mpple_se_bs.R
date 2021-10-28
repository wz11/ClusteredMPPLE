
library(survival)
library(data.table)
library(geepack)

clustered_mpple_se_bs <- function(data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, nboot = 1000, ...){

  b1 <- length(all.vars(formula21)) - 2
  b2 <- length(all.vars(formula22)) - 2
  x <- sort(unique(data$x))
  t <- sort(unique(data$x[data$c>0 & data$r==1]))
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
    data_temp <- rbindlist(lapply(cluster, function(x) data[data$clusterid == x,]))

    result_temp1 <- clustered_mpple_est(data, formula1 =  formula11, formula2 = formula21, cause = 1, w = TRUE)
    result_temp2 <- clustered_mpple_est(data, formula1 =  formula12, formula2 = formula22, cause = 2, w = TRUE)

    cif_temp1 <- clustered_cif_est(est1 = result_temp1, est2 = result_temp2, cause = 1, x)
    cif_temp2 <- clustered_cif_est(est1 = result_temp1, est2 = result_temp2, cause = 2, x)

    beta_boot_cause1[i,] <- result_temp1$beta
    beta_boot_cause2[i,] <- result_temp2$beta

    crp_boot_cause1[i,] <- result_temp1$crp
    crp_boot_cause2[i,] <- result_temp2$crp

    CIF_boot_cause1[i,] <- cif_temp1
    CIF_boot_cause2[i,] <- cif_temp2

  }

  out <- list(beta_boot_cause1 = beta_boot_cause1, beta_boot_cause2 = beta_boot_cause2,
              crp_boot_cause1 = crp_boot_cause1, crp_boot_cause2 = crp_boot_cause2,
              CIF_boot_cause1 = CIF_boot_cause1, CIF_boot_cause2 = CIF_boot_cause2)

  return(out)
}
