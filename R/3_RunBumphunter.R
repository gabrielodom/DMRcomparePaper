#' Return Results from the \code{bumphunter} Function
#'
#' @description A wrapper function for the Bumphunter method as implemented in
#'    the \code{bumphunter} package, called internally by the
#'    \code{\link{WriteBumphunterResults}} function.
#'
#' @param betaVals_mat A matrix of beta values returned in the first entry of
#'    the output from the \code{SimulateData} function, ordered by the CpGs. 
#'    Note this dataset inlcudes all CpGs on the array. 
#'    
#' @param labels_fct A factor vector of subject class labels. These should
#'    match the observations contained in the columns of the \code{betaVals_mat}
#'    matrix. Defaults to seven \code{"Tumor"} followed by seven \code{"Normal"}
#'    samples.
#'    
#' @param chromos_char A character vector for the chromosomes on which the CpGs
#'    are located
#' 
#' @param chromPosit_num A numeric vector for the locations of the CpGs
#' 
#' @param cpgLocation_df An annotation table that indicates locations of CpGs.
#'    This data frame has CpG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CpG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set.
#'    
#' @param pickCutoffQ_num The quantile used for picking the cutoff using the
#'    permutation distribution, passed to the \code{\link[bumphunter]{bumphunter}}
#'    function.
#'    
#' @param maxGap_int The maximum location gap, passed to the
#'    \code{\link[bumphunter]{bumphunter}} function. 
#'    
#' @param B_int An integer denoting the number of resamples to use when
#'    computing null distributions, passed to the
#'    \code{\link[bumphunter]{bumphunter}} function.
#'    
#' @param numCores The number of computing cores for parallel execution, passed
#'    to the \code{\link[doParallel]{registerDoParallel}} function. Defaults to
#'    one less than the number of cores available on your machine, as detected
#'    via the \code{\link[parallel]{detectCores}} function.
#'    
#' @param dmr.sig.threshold Regions with DMR p-value less than
#'    \code{dmr.sig.threshold} are selected for the output. 
#'    
#' @param min.cpgs Minimum number of CpGs. Regions with at least \code{min.cpgs}
#'    are selected for the output. Defaults to 5.
#'
#' @return A list of two elements: a data frame of \code{bumphunter} results
#'    and the computing time for the bumphunter method.
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
#'                             Aclusters_df = startEndCPG_df,
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
