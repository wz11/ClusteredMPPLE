

#' Obtain confidence band for cumulative incidence functions based on bootstrap variance
#'
#' @description This function provides equal precision bands and Hall-Wellner bands for cumulative incidence functions.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param CIF a list of time points and cumulative incidence functions.
#' @param CIF_boot a matrix of cumulative incidence function bootstrap samples.
#' @param range_lower lower limit for cumulative incidence functions.
#' @param range_upper upper limit for cumulative incidence functions.
#' @param nc number of clusters.
#' @param ... for future methods.
#'
#' @return A list with components:
#' \item{LB_EP}{a vector for lower limit of equal precision bands}
#' \item{UB_EP}{a vector for upper limit of equal precision bands}
#' \item{LB_HW}{a vector for lower limit of Hall-Wellner bands}
#' \item{UB_HW}{a vector for upper limit of Hall-Wellner bands}
#' \item{range}{a vector for range of confidence bands}
#'

clustered_cif_band_bs <- function(CIF, CIF_boot, range_lower, range_upper, nc, ...){

  x <- CIF$x
  range <- (x >= range_lower) & (x <= range_upper)

  cif_x <- CIF$CIF
  V_cif_x <- apply(CIF_boot, 2, stats::var)
  qEP_x <- cif_x * log(cif_x) / sqrt(V_cif_x)
  qHW_x <- cif_x * log(cif_x) / (1 + V_cif_x)

  W_t <- sqrt(nc)* (t(CIF_boot)-CIF$CIF)

  B_t_EP <- abs(qEP_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t[range,])
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca_EP <- as.numeric(stats::quantile(B_t_EP, 0.95, na.rm = TRUE))

  B_t_HW <- abs(qHW_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t[range,])
  B_t_HW <- apply(B_t_HW, 2, max, na.rm = TRUE)
  ca_HW <- as.numeric(stats::quantile(B_t_HW, 0.95, na.rm = TRUE))

  LB_EP <- exp(-exp(log(-log(cif_x)) - ca_EP/(sqrt(nc) * qEP_x)))
  UB_EP <- exp(-exp(log(-log(cif_x)) + ca_EP/(sqrt(nc) * qEP_x)))

  LB_HW <- exp(-exp(log(-log(cif_x)) - ca_HW/(sqrt(nc) * qHW_x)))
  UB_HW <- exp(-exp(log(-log(cif_x)) + ca_HW/(sqrt(nc) * qHW_x)))

  out <- list(LB_EP=LB_EP,UB_EP=UB_EP,LB_HW=LB_HW,UB_HW=UB_HW,range=range)
  return(out)
}
