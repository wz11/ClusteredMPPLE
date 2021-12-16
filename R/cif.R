

#' Obtain Cumulative Incidence Function and Confidence Band
#'
#' @description \code{cif} method for class \code{ccr_smreg}. It provides the cumulative incidence function with given covariates and provide equal precision bands and Hall-Wellner bands.
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param object an object of class \code{ccr_smreg}, which is a result of a call to \code{ccr_smreg}.
#' @param band logical value: if TRUE, it will calculate and return the 95\% simutanous confidence bands.
#' @param sims a number of simulations used for calculating 95\% simutanous confidence bands.
#' @param ... for future methods.
#'
#' @return An object of class \code{cif.ccr_smreg} with components:
#' \item{x}{time points}
#' \item{CIF1}{the cumulative incidence function for cause 1}
#' \item{CIF2}{the cumulative incidence function for cause 2}
#' \item{CIF1_band}{a list of 95\% simutanous confidence bands for cause 1}
#' \item{CIF2_band}{a list of 95\% simutanous confidence bands for cause 2}
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("simulated_data")
#' fit <- ccr_smreg(data=simulated_data, formula1 =  y ~ x + Z1 + Z2,
#' formula2 = Surv(x, d) ~ Z1 + Z2, ics.weight = TRUE, var.method = "BS", nboot = 100)
#' ## Continuing the ccr_smreg(...) example
#' pfit <- cif.ccr_smreg(object = fit, band = TRUE, sims = 1000)
#' plot(pfit)
#' }
#' @export

cif.ccr_smreg <- function(object, band = TRUE, sims = 1,...){

  var.method <- attr(object, "var.method")
  var.method <- match.arg(var.method, c("None", "CF", "BS"))

  cif1 <- clustered_cif_est(est1 = object$est1, est2 = object$est2, cause = 1, t = object$x)
  cif2 <- clustered_cif_est(est1 = object$est1, est2 = object$est2, cause = 2, t = object$x)

  if (band == TRUE){

    if(var.method == "None"){

      stop("x has attribute var.method = 'None'")

    } else if(var.method == "CF"){

      cif_band1 <- clustered_cif_band_cf(se1 = object$se1, se2 = object$se2, x = object$x, range_lower = object$Q1, range_upper = object$Q9, cause = 1, sims = sims)
      cif_band2 <- clustered_cif_band_cf(se1 = object$se1, se2 = object$se2, x = object$x, range_lower = object$Q1, range_upper = object$Q9, cause = 2, sims = sims)

    } else if(var.method == "BS"){

      cif_band1 <- clustered_cif_band_bs(CIF=cif1, CIF_boot=object$se1$CIF_boot, range_lower= object$Q1, range_upper= object$Q9, nc = object$nc)
      cif_band2 <- clustered_cif_band_bs(CIF=cif2, CIF_boot=object$se2$CIF_boot, range_lower= object$Q1, range_upper= object$Q9, nc = object$nc)

    }

    out <- list(x = object$x, CIF1 = cif1$CIF, CIF2 =  cif2$CIF, CIF1_band = cif_band1, CIF2_band = cif_band2)

  } else if (band == FALSE){
    out <- list(x = object$x, CIF1 = cif1$CIF, CIF2 =  cif2$CIF, CIF1_band = NULL, CIF2_band = NULL)
  }

  class(out) <- "cif.ccr_smreg"
  attr(out, "band") <- band
  return(out)
}


#' Plot a \code{cif.ccr_smreg} object
#'
#' @description Plot the cumulative incidence functions and 95\% simutanous confidence bands from a \code{cif.ccr_smreg} object.
#'
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param x the result of a call to 'cif.ccr_smreg'
#' @param ... for future methods
#'
#' @examples
#' \dontrun{
#' plot(pfit)
#' }
#' @export


plot.cif.ccr_smreg <- function(x,...){

  band <- attr(x, "band")
  graphics::par(mfrow = c(1, 2))

  t <- x$x
  CIF1 <- x$CIF1
  CIF2 <- x$CIF2

  if(band == TRUE){

    range <- x$CIF1_band$range

    graphics::plot(t[range], CIF1[range], type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 1")
    graphics::lines(t[range], x$CIF1_band$LB_EP[range], lty = 3, lwd = 2)
    graphics::lines(t[range], x$CIF1_band$UB_EP[range], lty = 3, lwd = 2)
    graphics::lines(t[range], x$CIF1_band$LB_HW[range], lty = 5, lwd = 2)
    graphics::lines(t[range], x$CIF1_band$UB_HW[range], lty = 5, lwd = 2)

    graphics::plot(t[range], CIF2[range], type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 2")
    graphics::lines(t[range], x$CIF2_band$LB_EP[range], lty = 3, lwd = 2)
    graphics::lines(t[range], x$CIF2_band$UB_EP[range], lty = 3, lwd = 2)
    graphics::lines(t[range], x$CIF2_band$LB_HW[range], lty = 5, lwd = 2)
    graphics::lines(t[range], x$CIF2_band$UB_HW[range], lty = 5, lwd = 2)

    }else if(band == FALSE){

      graphics::plot(t, CIF1, type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 1")
      graphics::plot(t, CIF2, type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 2")

  }
}
