
#' pseudo_HIV_data data set
#'
#' @description We generate sample HIV dataset for this package
#'
#' @docType data
#' @format A data.frame with 7 columns:
#' \describe{
#'   \item{clusterid}{Cluster ID, number from 1 to the number of clusters}
#'   \item{x}{Observed time}
#'   \item{c}{Cause of failure, 1 = cause 1, 2 = cause 2, 0 = right censoring, 99 = missing cause of failure}
#'   \item{Sex}{Covariate Sex}
#'   \item{Age}{Covariate Age}
#'   \item{CD4}{Covariate CD4}
#'   \item{HIV}{Covariate HIV}
#' }
#'
#' @source Simulated Data
#' @examples
#' library(ClusteredMPPLE)
#' data(pseudo_HIV_data)
"pseudo_HIV_data"
