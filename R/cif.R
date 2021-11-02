
cif.ccr_smreg <- function(object, band = TRUE, sims = 1,...){

  var <- attr(object, "var")
  var <- match.arg(var, c("None", "CF", "BS"))

  cif1 <- clustered_cif_est(est1 = object$est1, est2 = object$est2, cause = 1, t = object$x)
  cif2 <- clustered_cif_est(est1 = object$est1, est2 = object$est2, cause = 2, t = object$x)

  if (band == TRUE){

    if(var == "None"){

      stop("x has attribute var = 'None'")

    } else if(var == "CF"){

      cif_band1 <- clustered_cif_band_cf(se1 = object$se1, se2 = object$se2, x = object$x, range_lower = object$Q1, range_upper = object$Q9, cause = 1, sims = sims)
      cif_band2 <- clustered_cif_band_cf(se1 = object$se1, se2 = object$se2, x = object$x, range_lower = object$Q1, range_upper = object$Q9, cause = 2, sims = sims)

    } else if(var == "BS"){

      cif_band1 <- clustered_cif_band_bs(CIF=cif1, CIF_boot=object$se1$CIF_boot, range_lower= object$Q1, range_upper= object$Q9, nc = object$nc)
      cif_band2 <- clustered_cif_band_bs(CIF=cif2, CIF_boot=object$se1$CIF_boot, range_lower= object$Q1, range_upper= object$Q9, nc = object$nc)

    }

    out <- list(x = object$x, CIF1 = cif1$CIF, CIF2 =  cif2$CIF, CIF1_band = cif_band1, CIF2_band = cif_band2)

  } else if (band == FALSE){
    out <- list(x = object$x, CIF1 = cif1$CIF, CIF2 =  cif2$CIF)
  }

  class(out) <- "cif.ccr_smreg"
  attr(out, "band") <- band
  return(out)
}


plot.cif.ccr_smreg <- function(object,...){

  band <- attr(object, "band")
  par(mfrow = c(1, 2))

  x <- object$x
  CIF1 <- object$CIF1
  CIF2 <- object$CIF2

  if(band == TRUE){

    range <- object$CIF1_band$range

    graphics::plot(x[range], CIF1[range], type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 1",...)
    graphics::lines(x[range], object$CIF1_band$LB_EP[range], lty = 3, lwd = 2)
    graphics::lines(x[range], object$CIF1_band$UB_EP[range], lty = 3, lwd = 2)
    graphics::lines(x[range], object$CIF1_band$LB_HW[range], lty = 5, lwd = 2)
    graphics::lines(x[range], object$CIF1_band$UB_HW[range], lty = 5, lwd = 2)

    graphics::plot(x[range], CIF2[range], type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 2",...)
    graphics::lines(x[range], object$CIF2_band$LB_EP[range], lty = 3, lwd = 2)
    graphics::lines(x[range], object$CIF2_band$UB_EP[range], lty = 3, lwd = 2)
    graphics::lines(x[range], object$CIF2_band$LB_HW[range], lty = 5, lwd = 2)
    graphics::lines(x[range], object$CIF2_band$UB_HW[range], lty = 5, lwd = 2)

    }else if(band == FALSE){

      graphics::plot(x, CIF1, type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 1",...)
      graphics::plot(x, CIF2, type="l", ylim = c(0,1), lty = 1, lwd = 3, xlab = "Time", ylab = "Cumulative incidence function", main = "Cause 2",...)

  }
}
