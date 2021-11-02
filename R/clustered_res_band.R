

clustered_res_band <- function(crp, crp_boot, t, range_lower, range_upper, nc){

  range_r <- (t >= range_lower) & (t <= range_upper)
  res_x <- crp
  W_t <- sqrt(nc)* (t(crp_boot) - crp)

  B_t_EP <- abs(W_t[range_r,])
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca <- as.numeric(quantile(B_t_EP, 0.95, na.rm = TRUE)) / sqrt(nc)

  p.value <- mean(B_t_EP > max(abs(crp)))
  out <- list(ca = ca, p.value = p.value, range_r = range_r)
  return(out)

}
