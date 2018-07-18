#' Return Results from the \code{champ.DMR} Function
#'
#' @description A wrapper function for the ProbeLasso method, called internally
#'    by the \code{\link{WriteProbeLassoResults}} function. This function calls
#'    the \code{\link[ChAMP]{champ.DMR}} function to perform the ProbeLasso
#'    method calculations.
#'
#' @param betaVals_mat A matrix of beta values returned in the second entry of
#'    the output from the \code{SimulateData} function
#' @param labels_fct A factor vector of subject class labels. These should
#'    match the observations contained in the columns of the \code{betaVals_mat}
#'    matrix. Defaults to \code{factor(c(rep("Tumor", 7), rep("Normal", 7)))}
#' @param cpgLocation_df A data frame matching chromosomes to CPG names and
#'    locations. This is given in the \code{cpgLocation_df} data set.
#' @param adjPvalProbe_num The minimum threshold of significance for probes to
#'    be includede in DMRs, passed to the \code{\link[ChAMP]{champ.DMR}}
#'    function.
#' @param meanLassoRadius_int Radius around each DMP to detect DMR, passed to
#'    the \code{\link[ChAMP]{champ.DMR}} function.
#' @param minDmrSep_int The minimum seperation (bp) between neighbouring DMRs,
#'    passed to the \code{\link[ChAMP]{champ.DMR}} function.
#' @param nCores How many cores should be used to perform calculations? Defaults
#'    to 1. Note that this function should be called from within the
#'    \code{\link{WriteProbeLassoResults}} function, which is already written
#'    in parallel. If this function is executed directly (not from within this
#'    function), then this argument is passed to the \code{cores} argument of
#'    the \code{\link[ChAMP]{champ.DMR}} function.
#' @param dmr.sig.threshold Significance level to select regions (with
#'    \code{dmr.pval} less than the specified value) passed to the internal
#'    \code{\link{StandardizeOutput}} function
#' @param min.cpgs The minimum number of CPGs necessary to consider a result
#'    significant, passed to the internal \code{\link{StandardizeOutput}}
#'    function. Defaults to 5.
#'
#' @return A list of two elements: a data frame of \code{champ.DMR} results that
#'    have been standardized by the \code{\link{StandardizeOutput}} function and
#'    the computing time for the ProbeLasso method.
#'
#' @import ChAMPdata
#'
#' @importFrom ChAMP champ.DMR
#'
#' @export
#'
#' @examples
#' # Called internally by the WriteProbeLassoResults() function.
#' \dontrun{
#'    data("betaVals_mat")
#'    data("cpgLocation_df")
#'    data("startEndCPG_df")
#'
#'    treat_ls <- SimulateData(beta_mat = betaVals_mat,
#'                             AclustCPG_df = startEndCPG_df,
#'                             delta_num = 0.4,
#'                             seed_int = 100)
#'    class_fct <- factor(c(rep("Tumor", 7), rep("Normal", 7)))
#'
#'    RunProbeLasso(
#'      betaVals_mat = treat_ls$simBetaVals_df,
#'      labels_fct = class_fct,
#'      cpgLocation_df = cpgLocation_df,
#'      adjPvalProbe_num = 0.05,
#'      meanLassoRadius_int = 1000,
#'      minDmrSep_int = 1000
#'    )
#' }
RunProbeLasso <- function(betaVals_mat,
                          labels_fct = factor(c(rep("Tumor", 7),
                                                rep("Normal", 7))),
                          cpgLocation_df,
                          adjPvalProbe_num,
                          meanLassoRadius_int,
                          minDmrSep_int,
                          nCores = 1,
                          dmr.sig.threshold = 0.05,
                          min.cpgs = 5) {

  ###  Calculate ProbeLasso Results  ###
  ptm <- proc.time()

  # Requires ChAMPdata::probe.features
  probeLasso_out <- tryCatch(

    suppressMessages(

      champ.DMR(beta = as.matrix(betaVals_mat),
                pheno = labels_fct,
                method = "ProbeLasso",
                minProbes = 3,
                adjPvalDmr = 1,
                cores = nCores,
                meanLassoRadius = meanLassoRadius_int,
                minDmrSep = minDmrSep_int,
                adjPvalProbe = adjPvalProbe_num,
                PDFplot = FALSE,
                Rplot = FALSE)

    ),
    error = function(e1){ NULL }

  )

  elapsedtime <- proc.time() - ptm

  ###  Extract Results  ###
  # extract results if any predicted cluster is identified from ProbeLasso
  if(!is.null(probeLasso_out)){

    probeLassoOut_df <- probeLasso_out$ProbeLassoDMR

    results_df <- StandardizeOutput(
      methodOut_df = probeLassoOut_df,
      method = "ProbeLasso",
      cpgLocation_df = cpgLocation_df,
      dmr.sig.threshold = dmr.sig.threshold,
      min.cpgs = min.cpgs
    )

  } else {
    results_df <- NULL
  }

  ###  Return  ###
  list(results_df, elapsedtime[3])

}
