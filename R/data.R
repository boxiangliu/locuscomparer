#' Internal database connection configuration.
#'
#' A one-row data.frame holding the MySQL connection parameters (host, user,
#' password, port and database name) used internally by
#' \code{\link{get_position}} and \code{\link{retrieve_LD}} to query genomic
#' positions and pairwise LD from the LocusCompare database.
#'
#' @docType data
#' @keywords internal
#' @format A data.frame with 1 row and 5 columns.
#' @name config
"config"
