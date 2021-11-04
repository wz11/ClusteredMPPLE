
#' Obtain goodness-of-fit for cumulative residual process
#'
#' @description This function provides goodness-of-fit bands for cumulative residual process.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param crp a vector of cumulative residual process.
#' @param crp_boot a matrix of cumulative residual process bootstrap samples.
#' @param t a vector of time points.
#' @param range_lower lower limit for cumulative residual process.
#' @param range_upper upper limit for cumulative residual process.
#' @param nc number of clusters.
#' @param ... for future arguments.
#'
#' @return A list with components:
#' \item{ca}{a number for limit of goodness-of-fit bands}
#' \item{p.value}{a number for p-value of goodness-of-fit test}
#' \item{range_r}{a vector for range of goodness-of-fit bands}
#'

clustered_res_band <- function(crp, crp_boot, t, range_lower, range_upper, nc,...){

  range_r <- (t >= range_lower) & (t <= range_upper)
  res_x <- crp
  W_t <- sqrt(nc)* (t(crp_boot) - crp)

  B_t_EP <- abs(W_t[range_r,])
  B_t_EP <- apply(B_t_EP, 2, max, na.rm = TRUE)
  ca <- as.numeric(stats::quantile(B_t_EP, 0.95, na.rm = TRUE)) / sqrt(nc)

  p.value <- mean(B_t_EP > max(abs(crp)))
  out <- list(ca = ca, p.value = p.value, range_r = range_r)
  return(out)

}
