#' mypaired
#'
#' A JAGS script for Bayesian analysis for paired data
#'
#' @param Data data source
#'
#' @return JAGS stuff
#' @export
#'
#' @examples \dontrun{mypaird(df)}
mypaired <- function(Data) {
  ys = Data$STANDARD
  yn = Data$HUFFMAN

  Ntotal = length(yn)
  y <- structure(c(ys, yn), .Dim=c(Ntotal,2))

  dataList = list(y=y, Ntotal=Ntotal)
  source(system.file("jags/LS11/DBDA2E-utilities.R", package = "DSALS11"), local = TRUE)
  source(system.file("jags/LS11/Jags-PairedSampleScriptMV-wayne.R", package = "DSALS11"), local = TRUE)
}
