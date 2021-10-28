
clustered_cif_band_cf <- function(se1, se2, range_lower, range_upper, times = c(0.2, 0.4, 0.8, 1.6), nc = 50, cause = 1, sims = 1000){

  H <- list()
  H[[1]] <- se1$H
  H[[2]] <- se2$H
  phi <- list()
  phi[[1]] <- se1$phi
  phi[[2]] <- se2$phi
  dphi <- list()
  dphi[[1]] <- se1$dphi
  dphi[[2]] <- se2$dphi
  x <- se1$time
  range <- (x >= range_lower) & (x <= range_upper)

  S <- exp(-H[[1]] - H[[2]])
  S_minus <- c(1, S[1 : (length(S) - 1)])
  Hj <- H[[cause]]
  Hj_minus <- c(0, Hj[1 : (length(Hj) - 1)])

  cif_j <- function(t){
    sum(S_minus*(Hj - Hj_minus)*(x <= t))
  }

  CIF <- sapply(times, cif_j, simplify = TRUE)
  phi1_minus <- rbind(rep(0, times = dim(phi[[1]])[2]),
                      phi[[1]][1 : (dim(phi[[1]])[1] - 1), ])
  phi2_minus <- rbind(rep(0, times = dim(phi[[2]])[2]),
                      phi[[2]][1 : (dim(phi[[2]])[1] - 1), ])

  zi_t <- function(t){
    colSums(S_minus*dphi[[cause]]*(x <= t)) -
      colSums((S_minus*(Hj - Hj_minus))*(phi1_minus + phi2_minus)*(x <= t))
  }

  zi <- sapply(times, zi_t, simplify = TRUE)
  nc <- dim(zi)[1]
  V_CIF <- sqrt(colMeans(zi^2)/nc)

  zi <- t(sapply(x, zi_t, simplify = TRUE))
  cif_x <- sapply(x, cif_j, simplify = TRUE)
  V_cif_x <- rowMeans(zi^2)

  qEP_x <- cif_x * log(cif_x) / sqrt(V_cif_x)
  qHW_x <- cif_x * log(cif_x) / (1 + V_cif_x)

  xi <- matrix(rnorm(n = nc*sims), ncol = sims)
  W_t <- (zi[range, ] %*% xi) / sqrt(nc)

  B_t_EP <- abs(qEP_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca_EP <- as.numeric(quantile(B_t_EP, 0.95, na.rm = TRUE))

  B_t_HW <- abs(qHW_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
  B_t_HW <- apply(B_t_HW, 2, max, na.rm = TRUE)
  ca_HW <- as.numeric(quantile(B_t_HW, 0.95, na.rm = TRUE))

  LB_EP <- exp(-exp(log(-log(cif_x)) - ca_EP/(sqrt(nc) * qEP_x)))
  UB_EP <- exp(-exp(log(-log(cif_x)) + ca_EP/(sqrt(nc) * qEP_x)))


  LB_HW <- exp(-exp(log(-log(cif_x)) - ca_HW/(sqrt(nc) * qHW_x)))
  UB_HW <- exp(-exp(log(-log(cif_x)) + ca_HW/(sqrt(nc) * qHW_x)))

  out <- list(times = times, V_CIF = V_CIF, x = x, LB_EP = LB_EP, UB_EP = UB_EP, LB_HW = LB_HW, UB_HW = UB_HW)
  return(out)

}
