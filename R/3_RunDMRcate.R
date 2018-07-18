#' Return Results from the \code{dmrcate} Function
#'
#' @description A wrapper function for the DMRcate method from the
#'    \code{DMRcate} package, called internally by the
#'    \code{\link{WriteDMRcateResults}} function.
#'
#' @param betaVals_mat A matrix of beta values returned in the second entry of
#'    the output from the \code{SimulateData} function.
#' @param labels_fct A factor vector of subject class labels. These should
#'    match the observations contained in the columns of the \code{betaVals_mat}
#'    matrix. Defaults to \code{factor(c(rep("Tumor", 7), rep("Normal", 7)))}
#' @param cpgLocation_df A data frame matching chromosomes to CPG names and
#'    locations. This is given in the \code{cpgLocation_df} data set.
#' @param lambda_int Gaussian kernel bandwidth for smoothed-function estimation
#'    in the called \code{\link[DMRcate]{dmrcate}} function.
#' @param C_int Scaling factor for bandwidth in the internal call to the
#'    \code{\link[DMRcate]{dmrcate}} function
#' @param nCores How many cores should be used to perform calculations? Defaults
#'    to 1. Note that this function should be called from within the
#'    \code{\link{WriteDMRcateResults}} function, which is already written in
#'    parallel. Further note that the \code{DMRcate} package (as of version
#'    1.16.0), does not support parallelization in Windows environments.
#' @param dmr.sig.threshold Significance level to select regions (with
#'    \code{dmr.pval} less than the specified value) passed to the internal
#'    \code{\link{StandardizeOutput}} function.
#' @param min.cpgs The minimum number of CPGs necessary to consider a result
#'    significant, passed to the internal \code{\link{StandardizeOutput}}
#'    function. Defaults to 5.
#' @param genome Reference genome for annotating DMRs with promoter overlaps,
#'    passed to the \code{\link[DMRcate]{extractRanges}} function. Can be one
#'    of \code{"hg19"}, \code{"hg38"}, or \code{"mm10"}. Defaults to
#'    \code{"hg19"}.
#'
#'
#' @return A list of two elements: a data frame of \code{dmrcate} results
#'    that have been standardized by the \code{\link{StandardizeOutput}}
#'    function and the computing time for the DMRcate method.
#'
#' @import DMRcatedata
#'
#' @importFrom minfi logit2
#' @importFrom stats model.matrix
#' @importFrom DMRcate cpg.annotate
#' @importFrom DMRcate dmrcate
#' @importFrom DMRcate extractRanges
#'
#' @export
#'
#' @examples
#' # Called internally by the WriteDMRcateResults() function.
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
#'    RunDMRcate(
#'      betaVals_mat = treat_ls$simBetaVals_df,
#'      labels_fct = class_fct,
#'      cpgLocation_df = cpgLocation_df,
#'      lambda_int = 500, C_int = 5
#'    )
#' }
RunDMRcate <- function(betaVals_mat,
                       labels_fct = factor(c(rep("Tumor", 7),
                                             rep("Normal", 7))),
                       cpgLocation_df,
                       lambda_int, C_int,
                       nCores = 1,
                       dmr.sig.threshold = 0.05,
                       min.cpgs = 5, genome = "hg19"){


  ###  Calculate dmrcate Results  ###
  ptm <- proc.time()
  M_mat <- logit2(as.matrix(betaVals_mat))
  design_mat <- model.matrix(~labels_fct)

  # This takes 53.22585 sec
  myannotation <- cpg.annotate("array", M_mat, what = "M",
                               arraytype = "450K",
                               analysis.type = "differential",
                               design = design_mat,
                               coef = 2, fdr = 0.05)

  if(Sys.info()['sysname'] == "Windows"){
    nCores <- 1
    warning("DMRcate v 1.16.0 does not support parallelization on Windows. Executing serially.")
  }

  # This takes 8.167 sec: 13% of computing time.
  dmrcate_out <- tryCatch(

    # Parallel no supported on Windows
    dmrcate(myannotation, lambda = lambda_int, C = C_int, mc.cores = nCores),
    error = function(e1){ NULL }

  )

  elapsedtime <- proc.time() - ptm


  ###  Extract Results  ###

  # If any predicted cluster is identified from DMRcate
  if(!is.null(dmrcate_out)){

    # Requires DMRcatedata::dmrcatedata
    dmrcateOut_df <- data.frame(extractRanges(dmrcate_out, genome = genome))

    results_df <- StandardizeOutput(
      methodOut_df = dmrcateOut_df,
      method = "DMRcate",
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
