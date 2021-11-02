

crp.ccr_smreg <- function(object,...){
  var <- attr(object, "var")

  if (var != "BS"){
    stop("x does not have attribute var = 'BS'")
  }

  res1 <- clustered_res_band(crp = object$est1$crp, crp_boot = object$se1$crp_boot, t = object$t, range_lower = object$q1, range_upper = object$q9, nc = object$nc)
  res2 <- clustered_res_band(crp = object$est2$crp, crp_boot = object$se2$crp_boot, t = object$t, range_lower = object$q1, range_upper = object$q9, nc = object$nc)

  out <- list(res1_band = res1, res2_band = res2, t = object$t, res1 = object$est1$crp, res2 = object$est2$crp)
  class(out) <- "crp.ccr_smreg"
  return(out)
}


plot.crp.ccr_smreg <- function(object, cause,...){
  ca1 <- object$res1_band$ca
  ca2 <- object$res2_band$ca
  pval1 <- object$res1_band$p.value
  pval2 <- object$res2_band$p.value
  range_r <- object$res1_band$range_r
  t <- object$t

  xx <- c(t[range_r], rev(t[range_r]))
  yy1 <- c(rep(ca1*100,sum(range_r)), rep(-ca1*100,sum(range_r)))
  yy2 <- c(rep(ca2*100,sum(range_r)), rep(-ca2*100,sum(range_r)))
  par(mfrow=c(1,1))

  if (cause == 1){

    graphics::plot(t[range_r], object$res1[range_r]*100, type = "l", ylim = c(-ca1,ca1)*120, lty = 1, lwd = 2, xlab = "Time", ylab = expression("Residual Process ("%*%"100)"), main = "Cause 1",...)
    graphics::abline(h=0,...)
    graphics::polygon(xx, yy2, col = adjustcolor("gray",alpha.f=0.4), border = adjustcolor("gray",alpha.f=0.4),...)
    graphics::text(max(t[range_r])*0.8, -ca1*100, paste("p-value=",pval2), cex = 1,...)

  } else if (cause == 2){

    graphics::plot(t[range_r], object$res2[range_r]*100, type = "l", ylim = c(-ca2,ca2)*120, lty = 1, lwd = 2, xlab = "Time", ylab = expression("Residual Process ("%*%"100)"), main = "Cause 2",...)
    graphics::abline(h=0,...)
    graphics::polygon(xx, yy2, col = adjustcolor("gray",alpha.f=0.4), border = adjustcolor("gray",alpha.f=0.4),...)
    graphics::text(max(t[range_r])*0.8, -ca2*100, paste("p-value=", pval2), cex = 1, ...)

  }

}
