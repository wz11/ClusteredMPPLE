

#' Obtain confidence band for cumulative incidence functions based on close form variance
#'
#' @description This function provides equal precision bands and Hall-Wellner bands for cumulative incidence functions.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param se1 a list for the standard errors for cause 1.
#' @param se2 a list for the standard errors for cause 2.
#' @param x time points.
#' @param range_lower lower limit for cumulative incidence functions.
#' @param range_upper upper limit for cumulative incidence functions.
#' @param cause 1 or 2, the cause of failure to be evaluated.
#' @param sims a number of simulations used for calculating 95\% simutanous confidence bands.
#' @param ... for future methods.
#'
#' @return A list with components:
#' \item{V_CIF}{a vector of pointwise confidence intervals}
#' \item{LB_EP}{a vector for lower limit of equal precision bands}
#' \item{UB_EP}{a vector for upper limit of equal precision bands}
#' \item{LB_HW}{a vector for lower limit of Hall-Wellner bands}
#' \item{UB_HW}{a vector for upper limit of Hall-Wellner bands}
#' \item{range}{a vector for range of confidence bands}
#'

clustered_cif_band_cf <- function(se1, se2, x, range_lower, range_upper, cause = 1, sims = 1000,...){

  H <- list()
  H[[1]] <- se1$H
  H[[2]] <- se2$H
  phi <- list()
  phi[[1]] <- se1$phi
  phi[[2]] <- se2$phi
  dphi <- list()
  dphi[[1]] <- se1$dphi
  dphi[[2]] <- se2$dphi
  range <- (x >= range_lower) & (x <= range_upper)
  time <- se1$time
  S <- exp(-H[[1]] - H[[2]])
  S_minus <- c(1, S[1 : (length(S) - 1)])
  Hj <- H[[cause]]
  Hj_minus <- c(0, Hj[1 : (length(Hj) - 1)])

  cif_j <- function(t){
    sum(S_minus*(Hj - Hj_minus)*(time <= t))
  }


  phi1_minus <- rbind(rep(0, times = dim(phi[[1]])[2]),
                      phi[[1]][1 : (dim(phi[[1]])[1] - 1), ])
  phi2_minus <- rbind(rep(0, times = dim(phi[[2]])[2]),
                      phi[[2]][1 : (dim(phi[[2]])[1] - 1), ])

  zi_t <- function(t){
    colSums(S_minus*dphi[[cause]]*(time <= t)) -
      colSums((S_minus*(Hj - Hj_minus))*(phi1_minus + phi2_minus)*(time <= t))
  }

  zi <- sapply(x, zi_t, simplify = TRUE)
  nc <- dim(zi)[1]
  V_CIF <- sqrt(colMeans(zi^2)/nc)

  zi <- t(zi)
  cif_x <- sapply(x, cif_j, simplify = TRUE)
  V_cif_x <- rowMeans(zi^2)

  qEP_x <- cif_x * log(cif_x) / sqrt(V_cif_x)
  qHW_x <- cif_x * log(cif_x) / (1 + V_cif_x)

  xi <- matrix(stats::rnorm(n = nc*sims), ncol = sims)
  W_t <- (zi[range, ] %*% xi) / sqrt(nc)

  B_t_EP <- abs(qEP_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca_EP <- as.numeric(stats::quantile(B_t_EP, 0.95, na.rm = TRUE))

  B_t_HW <- abs(qHW_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t)
  B_t_HW <- apply(B_t_HW, 2, max, na.rm = TRUE)
  ca_HW <- as.numeric(stats::quantile(B_t_HW, 0.95, na.rm = TRUE))

  LB_EP <- exp(-exp(log(-log(cif_x)) - ca_EP/(sqrt(nc) * qEP_x)))
  UB_EP <- exp(-exp(log(-log(cif_x)) + ca_EP/(sqrt(nc) * qEP_x)))


  LB_HW <- exp(-exp(log(-log(cif_x)) - ca_HW/(sqrt(nc) * qHW_x)))
  UB_HW <- exp(-exp(log(-log(cif_x)) + ca_HW/(sqrt(nc) * qHW_x)))

  out <- list(V_CIF = V_CIF, LB_EP = LB_EP, UB_EP = UB_EP, LB_HW = LB_HW, UB_HW = UB_HW, range = range)
  return(out)

}
