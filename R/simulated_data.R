
#' simulated_data data set
#'
#' @description We simulated a sample dataset with 50 cluster size for this package
#'
#' @docType data
#' @format A data.frame with 5 columns:
#' \describe{
#'   \item{clusterid}{Cluster ID, number from 1 to the number of clusters}
#'   \item{x}{Observed time}
#'   \item{c}{Cause of failure, 1 = cause 1, 2 = cause 2, 0 = right censoring, 99 = missing cause of failure}
#'   \item{Z1}{Covariate Z1}
#'   \item{Z2}{Covariate Z2}
#' }
#'
#' @source Simulated Data
#' @examples
#' library(ClusteredMPPLE)
#' data(simulated_data)
"simulated_data"
