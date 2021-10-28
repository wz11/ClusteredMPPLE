
naive_cif_band_cf <- function(scen = scen, mpple_naive_temp1 = result_temp1_naive, mpple_naive_temp2 = result_temp2_naive, times=c(0.2, 0.4, 0.8, 1.6), cause = 1, band = 1, sims = 1000){
  H <- list()
  H[[1]] <- mpple_naive_temp1$H
  H[[2]] <- mpple_naive_temp2$H
  phi <- list()
  phi[[1]] <- mpple_naive_temp1$phi
  phi[[2]] <- mpple_naive_temp2$phi
  dphi <- list()
  dphi[[1]] <- mpple_naive_temp1$dphi
  dphi[[2]] <- mpple_naive_temp2$dphi
  x <- mpple_naive_temp1$x
  range <- mpple_naive_temp1$range

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
  n <- dim(zi)[1]
  V_CIF <- sqrt(colMeans(zi^2)/n)

  if(band == 1 & cause == 1){
    alpha <- 0.5
    zi <- t(sapply(x, zi_t, simplify = TRUE))
    cif_x <- sapply(x, cif_j, simplify = TRUE)
    V_cif_x <- rowMeans(zi^2)

    qEP_x <- cif_x * log(cif_x) / sqrt(V_cif_x)
    qHW_x <- cif_x * log(cif_x) / (1 + V_cif_x)

    xi <- matrix(rnorm(n = n*sims), ncol = sims)
    W_t <- (zi[range, ] %*% xi) / sqrt(n)

    B_t_EP <- abs(qEP_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
    B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
    ca_EP <- as.numeric(quantile(B_t_EP, 0.95, na.rm = TRUE))

    B_t_HW <- abs(qHW_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
    B_t_HW <- apply(B_t_HW, 2, max, na.rm = TRUE)
    ca_HW <- as.numeric(quantile(B_t_HW, 0.95, na.rm = TRUE))


    LB_EP <- exp(-exp(log(-log(cif_x)) - ca_EP/(sqrt(n) * qEP_x)))
    UB_EP <- exp(-exp(log(-log(cif_x)) + ca_EP/(sqrt(n) * qEP_x)))


    LB_HW <- exp(-exp(log(-log(cif_x)) - ca_HW/(sqrt(n) * qHW_x)))
    UB_HW <- exp(-exp(log(-log(cif_x)) + ca_HW/(sqrt(n) * qHW_x)))

    out <- list(LB_EP=LB_EP,UB_EP=UB_EP,LB_HW=LB_HW,UB_HW=UB_HW)
    return(out)
}
