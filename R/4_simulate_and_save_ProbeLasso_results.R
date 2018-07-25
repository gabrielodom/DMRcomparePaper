#' Calculate and Save ProbeLasso Method Results for Specified Design Points
#'
#' @description Given a set of design points, simulate appropriate DMR data and
#'    apply the ProbeLasso method (via the \code{\link[ChAMP]{champ.DMR}}
#'    function) to them (with parameters also within the design). Write the
#'    results to a file.
#'
#' @param beta_mat A beta value matrix for selected methylation samples from a
#'    450k methylation array with CPG IDs as the row names and sample IDs as
#'    the column names. An example is given in the \code{betaVals_mat} data set.
#' @param CPGs_df An annotation table that indicates locations of CpGs.
#'    This data frame has CPG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set.
#' @param Aclusters_df A data frame of beta values and CpG information for
#'    clusters of CpGs over a 450k methylation array. The rows correspond to the
#'    CPGs. The columns have information on the cluster number, chromosome,
#'    cluster start and end locations, and the beta values for each subject
#'    grouped by some clinical indicator (e.g. case v. control). An example is
#'    given in the \code{startEndCPG_df} data set.
#' @param parallel Should computing be completed over multiple computing cores?
#'    Defaults to \code{TRUE}.
#' @param numCores If \code{parallel}, how many cores should be used? Defaults
#'    to two less than the number of available cores (as calculated by the
#'    \code{\link[parallel]{detectCores}} function).
#' @param deltas_num A vector of treatment sizes: non-negative real numbers to
#'    add to the beta values within randomly-selected clusters for a single
#'    class of subjects. This artifically creates differentially-methylated
#'    regions (DMRs).
#' @param seeds_int A vector of seed values passed to the
#'    \code{\link[base]{Random}} function to enable reproducible results
#' @param pVals_num A vector of the minimum thresholds of significance for
#'    probes to be includede in DMRs, passed through the
#'    \code{\link{RunProbeLasso}} function to the \code{\link[ChAMP]{champ.DMR}}
#'    function.
#' @param aveLassoRad_int A vector of radii around each differential
#'    methylation position to detect DMR, passed to the
#'    \code{\link[ChAMP]{champ.DMR}} function.
#' @param minDmrSep_int A vector of the minimum seperation (bp) values between
#'    neighbouring DMRs, passed to the \code{\link[ChAMP]{champ.DMR}} function.
#' @param resultsDir Where should the results be saved? Defaults to
#'    \code{DMRcate_compare/}.
#' @param verbose Should the function print progress messages? Defaults to
#'    \code{TRUE} only if \code{parallel = FALSE}.
#'
#' @return Nothing. Saves output to a file in the specified results directory.
#'
#' @details This function creates matrices of all combinations of design points
#'    and all combinations of parameters. For each combination, this function
#'    executes the internal \code{\link{RunProbeLasso}} function and saves the
#'    results as a compressed \code{.RDS} file.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    data("betaVals_mat")
#'    data("cpgLocation_df")
#'    data("startEndCPG_df")
#'
#'    WriteProbeLassoResults(
#'      beta_mat = betaVals_mat,
#'      CPGs_df = cpgLocation_df,
#'      Aclusters_df = startEndCPG_df
#'    )
#' }
WriteProbeLassoResults <- function(beta_mat,
                                   CPGs_df,
                                   Aclusters_df,
                                   parallel = TRUE,
                                   numCores = detectCores() - 2,
                                   deltas_num = c(0, 0.025, 0.05, 0.10,
                                                  0.15, 0.20, 0.30, 0.40),
                                   seeds_int = c(100, 210, 330, 450, 680),
                                   pVals_num = c(0.001, 0.01, 0.05, 0.1),
                                   aveLassoRad_int = c(375, 700, 1000),
                                   minDmrSep_int = c(200, 250, 500, 750, 1000),
                                   resultsDir = "DMRcate_compare/",
                                   verbose = !parallel){

  dir.create(paste0("./", resultsDir), showWarnings = FALSE)

  ###  Data Simulation Outer Loop  ###
  designPts_mat <- expand.grid(deltas_num, seeds_int)
  paramsGrid_mat <- expand.grid(pVals_num, aveLassoRad_int, minDmrSep_int)

  for(i in 1:nrow(designPts_mat)){

    ###  Generate Data Set  ###
    delta <- designPts_mat[i, 1]
    seed  <- designPts_mat[i, 2]

    treatment_ls <- SimulateData(beta_mat = beta_mat,
                                 AclustCPG_df = Aclusters_df,
                                 delta_num = delta,
                                 seed_int = seed)
    betas_df <- treatment_ls$simBetaVals_df

    ###  Parallel Setup  ###
    if(!parallel){
      numCores <- 1
    }
    clust <- makeCluster(numCores)
    registerDoParallel(clust)


    ###  Inner Parameter Grid Search  ###
    foreach(j = 1:nrow(paramsGrid_mat),
            .packages = c("DMRcompare", "data.table")) %dopar% {

              ###  Calculate Method Output  ###
              adjPval   <- paramsGrid_mat[j, 1]
              mLassoRad <- paramsGrid_mat[j, 2]
              minDmrSep <- paramsGrid_mat[j, 3]

              res_ls <- RunProbeLasso(betaVals_mat = betas_df,
                                      cpgLocation_df = CPGs_df,
                                      adjPvalProbe_num = adjPval,
                                      meanLassoRadius_int = mLassoRad,
                                      minDmrSep_int = minDmrSep)

              ###  Define NULL Data  ###
              if(is.null(res_ls[[1]])){
                res_ls[[1]] <- data.frame(
                  dmr.chr      = NA,
                  dmr.start    = NA_integer_,
                  dmr.end      = NA_integer_,
                  seqnames     = NA,
                  start        = NA_integer_,
                  end          = NA_integer_,
                  width        = NA_integer_,
                  strand       = NA,
                  dmrNo        = NA_integer_,
                  dmrP         = NA_real_,
                  dmrpRank     = NA_integer_,
                  dmrChrom     = NA_character_,
                  dmrStart     = NA_integer_,
                  dmrEnd       = NA_integer_,
                  dmrSize      = NA_integer_,
                  dmrCoreStart = NA_integer_,
                  dmrCoreEnd   = NA_integer_,
                  dmrCoreSize  = NA_integer_,
                  ensemblID    = NA_character_,
                  geneSymbol   = NA_character_,
                  betaAv_Normal= NA_real_,
                  betaAV_Tumor = NA_real_,
                  dmr.pval     = NA_real_,
                  dmr.n.cpgs   = NA_integer_
                )
              }

              ###  Save Results  ###
              file_char <- paste0(
                resultsDir, "ProbeLassoResults_delta", delta,
                "_seed", seed,
                "_adjPvalProbe", adjPval,
                "_meanLassoRd", mLassoRad,
                "_minDmrSep", minDmrSep,
                ".RDS"
              )

              if(verbose){
                message("Saving results to file ", file_char, "\n")
              }

              saveRDS(res_ls, file = file_char)

            } # END for(j)

    stopCluster(clust)

  } # END for(i)

}

