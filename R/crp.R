


#' Obtain Cumulative Residual Process and Goodness-of-fit Confidence Band
#'
#' @description \code{crp} method for class \code{ccr_smreg}. It provides the cumulative residual process and goodness-of-fit bands.
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param object an object of class \code{ccr_smreg}, which is a result of a call to \code{ccr_smreg}.
#' @param ... for future methods.
#'
#' @return An object of class \code{crp.ccr_smreg} with components:
#' \item{t}{time points}
#' \item{res1}{the cumulative residual process for cause 1}
#' \item{res2}{the cumulative residual process for cause 1}
#' \item{res1_band}{a list of goodness-of-fit bands and p-value for cause 1}
#' \item{res2_band}{a list of goodness-of-fit bands and p-value for cause 2}
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("simulated_data")
#' fit <- ccr_smreg(data=simulated_data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2,
#' formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2,
#' w = TRUE, var.method = "BS", nboot = 100)
#' ## Continuing the ccr_smreg(...) example
#' rfit <- crp.ccr_smreg(object = fit)
#' plot(rfit)
#' }
#' @export
crp.ccr_smreg <- function(object,...){
  var <- attr(object, "var")

  if (var != "BS"){
    stop("x does not have attribute var = 'BS'")
  }

  res1_band <- clustered_res_band(crp = object$est1$crp, crp_boot = object$se1$crp_boot, t = object$t, range_lower = object$q1, range_upper = object$q9, nc = object$nc)
  res2_band <- clustered_res_band(crp = object$est2$crp, crp_boot = object$se2$crp_boot, t = object$t, range_lower = object$q1, range_upper = object$q9, nc = object$nc)

  out <- list(t = object$t, res1 = object$est1$crp, res2 = object$est2$crp, res1_band = res1_band, res2_band = res2_band)
  class(out) <- "crp.ccr_smreg"
  return(out)
}



#' Plot a \code{crp.ccr_smreg} object
#'
#' @description Plot the cumulative residual process from a \code{crp.ccr_smreg} object.
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param x the result of a call to 'crp.ccr_smreg'
#' @param cause 1 or 2, the cause of failure to be evaluated
#' @param ... for future methods
#' @examples
#' \dontrun{
#' plot(rfit)
#' }
#' @export

plot.crp.ccr_smreg <- function(x, cause = 1,...){
  ca1 <- x$res1_band$ca
  ca2 <- x$res2_band$ca
  pval1 <- x$res1_band$p.value
  pval2 <- x$res2_band$p.value
  range_r <- x$res1_band$range_r
  t <- x$t

  xx <- c(t[range_r], rev(t[range_r]))
  yy1 <- c(rep(ca1*100,sum(range_r)), rep(-ca1*100,sum(range_r)))
  yy2 <- c(rep(ca2*100,sum(range_r)), rep(-ca2*100,sum(range_r)))
  graphics::par(mfrow=c(1,1))

  if (cause == 1){

    graphics::plot(t[range_r], x$res1[range_r]*100, type = "l", ylim = c(-ca1,ca1)*120, lty = 1, lwd = 2, xlab = "Time", ylab = expression("Residual Process ("%*%"100)"), main = "Cause 1")
    graphics::abline(h=0)
    graphics::polygon(xx, yy2, col = grDevices::adjustcolor("gray",alpha.f=0.4), border = grDevices::adjustcolor("gray",alpha.f=0.4))
    graphics::text(max(t[range_r])*0.8, -ca1*100, paste("p-value=",pval1), cex = 1)

  } else if (cause == 2){

    graphics::plot(t[range_r], x$res2[range_r]*100, type = "l", ylim = c(-ca2,ca2)*120, lty = 1, lwd = 2, xlab = "Time", ylab = expression("Residual Process ("%*%"100)"), main = "Cause 2")
    graphics::abline(h=0)
    graphics::polygon(xx, yy2, col = grDevices::adjustcolor("gray",alpha.f=0.4), border = grDevices::adjustcolor("gray",alpha.f=0.4))
    graphics::text(max(t[range_r])*0.8, -ca2*100, paste("p-value=", pval2), cex = 1)

  }

}
