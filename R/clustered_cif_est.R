
#' Obtain Cumulative Incidence Function Estimates
#'
#' @description It provides the cumulative incidence function with given covariates for specific cause.
#' @author Wenxian Zhou, \email{wz11 at iu dot edu}
#' @author Giorgos Bakoyannis, \email{gbakogia at iu dot edu}
#'
#' @param est1 a list for the estimates for cause 1.
#' @param est2 a list for the estimates for cause 2.
#' @param cause 1 or 2, the cause of failure to be evaluated.
#' @param t time points.
#' @param ... for furthur arguments.
#'
#' @return Alist with components:
#' \item{x}{a vector for time points where cumulative incidence functions are evaluated}
#' \item{CIF}{a vector for the cumulative incidence function evaluated at time points x}
#'

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
  CIF_function <- stats::stepfun(time, c(0, CIF_l))

  CIF <- CIF_function(t)
  out <- list(x = t, CIF = CIF)
  return(out)
}
