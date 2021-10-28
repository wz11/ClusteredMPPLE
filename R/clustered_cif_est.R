


clustered_cif_est <- function(est1, est2, cause, t, ...){

  H <- list()
  H[[1]] <- est1$H
  H[[2]] <- est2$H
  time <- est1$time

  S <- exp(-H[[1]] - H[[2]])
  S_minus <- c(1, S[1 : (length(S) - 1)])

  Hl <- H[[cause]]
  Hl_minus <- c(0, Hl[1 : (length(Hl) - 1)])

  CIF_l <- cumsum(S_minus*(Hl - Hl_minus))
  CIF_function <- stepfun(time, c(0, CIF_l))

  CIF <- CIF_function(t)
  return(CIF)
}
