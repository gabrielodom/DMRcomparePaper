#' Return Results from the \code{bumphunter} Function
#'
#' @description A wrapper function for the Bumphunter method as implemented in
#'    the \code{bumphunter} package, called internally by the
#'    \code{\link{WriteBumphunterResults}} function.
#'
#' @param betaVals_mat A matrix of beta values returned in the second entry of
#'    the output from the \code{SimulateData} function, ordered by the CPGs.
#' @param labels_fct A factor vector of subject class labels. These should
#'    match the observations contained in the columns of the \code{betaVals_mat}
#'    matrix. Defaults to \code{factor(c(rep("Tumor", 7), rep("Normal", 7)))}
#' @param chromos_char A character vector with the chromosomes of each location
#' @param chromPosit_num A numeric vector representing the chromosomal position
#' @param cpgLocation_df An annotation table that indicates locations of CpGs.
#'    This data frame has CPG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set.
#' @param pickCutoffQ_num The quantile used for picking the cutoff using the
#'    permutation distribution, passed to the \code{\link[bumphunter]{bumphunter}}
#'    function.
#' @param maxGap_int The maximum location gap, passed to the
#'    \code{\link[bumphunter]{bumphunter}} function. This will be used to define
#'    the clusters of locations that are to be analyzed together via the
#'    \code{\link[bumphunter]{clusterMaker}} function.
#' @param B_int An integer denoting the number of resamples to use when
#'    computing null distributions, passed to the
#'    \code{\link[bumphunter]{bumphunter}} function.
#' @param numCores The number of computing cores for parallel execution, passed
#'    to the \code{\link[doParallel]{registerDoParallel}} function. Defaults to
#'    one less than the number of cores available on your machine, as detected
#'    via the \code{\link[parallel]{detectCores}} function.
#' @param dmr.sig.threshold Significance level to select regions (with
#'    \code{dmr.pval} less than the specified value) passed to the internal
#'    \code{\link{StandardizeOutput}} function.
#' @param min.cpgs The minimum number of CPGs before we consider a result
#'    significant, passed to the internal \code{\link{StandardizeOutput}}
#'    function. Defaults to 5.
#'
#' @return A list of two elements: a data frame of \code{bumphunter} results
#'    that have been standardized by the \code{\link{StandardizeOutput}}
#'    function and the computing time for the bumphunter method.
#'
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom bumphunter bumphunter
#'
#' @export
#'
#' @examples
#' # Called internally by the WriteBumphunterResults() function.
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
#'    RunBumphunter(
#'      betaVals_mat = treat_ls$simBetaVals_df,
#'      labels_fct = class_fct,
#'      cpgLocation_df = cpgLocation_df,
#'      pickCutoffQ_num = 0.95,
#'      maxGap_int = 250
#'    )
#' }
#'
RunBumphunter <- function(betaVals_mat,
                          labels_fct = factor(c(rep("Tumor", 7),
                                                rep("Normal", 7))),
                          chromos_char,
                          chromPosit_num,
                          cpgLocation_df,
                          pickCutoffQ_num,
                          maxGap_int,
                          B_int = 10,
                          numCores = detectCores() - 1,
                          dmr.sig.threshold = 0.05,
                          min.cpgs = 5) {

  ###  Calculate Bumphunter Results  ###
  ptm <- proc.time()
  M_mat <- logit2(as.matrix(betaVals_mat))
  design_mat <- model.matrix(~labels_fct)

  registerDoParallel(cores = numCores)
  bumphunt_out <- tryCatch(

    bumphunter(M_mat, design = design_mat,
               chr = chromos_char,
               pos = chromPosit_num,
               pickCutoff = TRUE,
               pickCutoffQ = pickCutoffQ_num,
               maxGap = maxGap_int,
               nullMethod = "permutation",
               B = B_int, verbose = TRUE,
               type = "M"),
    error = function(e1){ NULL }

  )

  # The bumphunter help documentation mentions that this function is necessary
  #   for  "cleanup, on Windows". I will have to test this package on Linux or
  #   Mac.
  if(Sys.info()[["sysname"]] == "Windows"){
    bumphunter:::foreachCleanup()
  }

  elapsedtime <- proc.time() - ptm


  ###  Extract Results  ###
  if(!is.null(bumphunt_out)){

    bumphuntOut_df <- bumphunt_out$table

    results_df <- StandardizeOutput(
      methodOut_df = bumphuntOut_df,
      method = "Bumphunter",
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
