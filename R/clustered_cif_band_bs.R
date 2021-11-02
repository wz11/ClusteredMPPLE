
clustered_cif_band_bs <- function(CIF, CIF_boot, range_lower, range_upper, nc){

  x <- CIF$x
  range <- (x >= range_lower) & (x <= range_upper)

  cif_x <- CIF$CIF
  V_cif_x <- apply(CIF_boot, 2, stats::var)
  qEP_x <- cif_x * log(cif_x) / sqrt(V_cif_x)
  qHW_x <- cif_x * log(cif_x) / (1 + V_cif_x)

  W_t <- sqrt(nc)* (t(CIF_boot)-CIF$CIF)

  B_t_EP <- abs(qEP_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t[range,])
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca_EP <- as.numeric(quantile(B_t_EP, 0.95, na.rm = TRUE))

  B_t_HW <- abs(qHW_x[range] / (log(cif_x[range]) * cif_x[range]) * W_t[range,])
  B_t_HW <- apply(B_t_HW, 2, max, na.rm = TRUE)
  ca_HW <- as.numeric(quantile(B_t_HW, 0.95, na.rm = TRUE))

  LB_EP <- exp(-exp(log(-log(cif_x)) - ca_EP/(sqrt(nc) * qEP_x)))
  UB_EP <- exp(-exp(log(-log(cif_x)) + ca_EP/(sqrt(nc) * qEP_x)))

  LB_HW <- exp(-exp(log(-log(cif_x)) - ca_HW/(sqrt(nc) * qHW_x)))
  UB_HW <- exp(-exp(log(-log(cif_x)) + ca_HW/(sqrt(nc) * qHW_x)))

  out <- list(LB_EP=LB_EP,UB_EP=UB_EP,LB_HW=LB_HW,UB_HW=UB_HW,range=range)
  return(out)
}
