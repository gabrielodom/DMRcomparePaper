#' Calculate and Save Bumphunter Method Results for Specified Design Points
#'
#' @description Given a set of design points, simulate appropriate DMR data and
#'    apply the \code{bumphunter} method to them (with parameters also within
#'    the design). Write the results to a file.
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
#'    \code{\link[parallel]{detectCores}} function). These cores are used
#'    internally by the \code{\link[bumphunter]{bumphunter}} function.
#' @param deltas_num A vector of treatment sizes: non-negative real numbers to
#'    add to the beta values within randomly-selected clusters for a single
#'    class of subjects. This artifically creates differentially-methylated
#'    regions (DMRs).
#' @param seeds_int A vector of seed values passed to the
#'    \code{\link[base]{Random}} function to enable reproducible results
#' @param cutoffQ_num A vector of quantiles used for picking the cutoff using
#'    the permutation distribution, passed through the call to the internal
#'    \code{\link{RunBumphunter}} call to \code{\link[bumphunter]{bumphunter}}.
#' @param maxGap_int A vector of maximum location gaps, passed to the
#'    \code{\link[bumphunter]{bumphunter}} function. These will be used to
#'    define the clusters of locations that are to be analyzed together via the
#'    \code{\link[bumphunter]{clusterMaker}} function.
#' @param resultsDir Where should the results be saved? Defaults to
#'    \code{DMRcate_compare/}.
#' @param verbose Should the function print progress messages? Defaults to
#'    \code{TRUE}.
#'
#' @return Nothing. Saves output to a file in the specified results directory.
#'
#' @details This function creates matrices of all combinations of design points
#'    and all combinations of parameters. For each combination, this function
#'    executes the internal \code{\link{RunBumphunter}} function and saves the
#'    results as a compressed \code{.RDS} file.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    data("betaVals_mat")
#'    data("cpgLocation_df")
#'    data("startEndCPG_df")
#'
#'    WriteBumphunterResults(
#'      beta_mat = betaVals_mat,
#'      CPGs_df = cpgLocation_df,
#'      Aclusters_df = startEndCPG_df
#'    )
#' }
WriteBumphunterResults <- function(beta_mat,
                                   CPGs_df,
                                   Aclusters_df,
                                   parallel = TRUE,
                                   numCores = detectCores() - 2,
                                   deltas_num = c(0, 0.025, 0.05, 0.10,
                                                  0.15, 0.20, 0.30, 0.40),
                                   seeds_int = c(100, 210, 330, 450, 680),
                                   cutoffQ_num = c(0.9, 0.95, 0.99),
                                   maxGap_int = c(200, 250, 500, 750, 1000),
                                   resultsDir = "DMRcate_compare/",
                                   verbose = TRUE){

  dir.create(paste0("./", resultsDir), showWarnings = FALSE)

  ###  Data Simulation Outer Loop  ###
  designPts_mat <- expand.grid(deltas_num, seeds_int)
  paramsGrid_mat <- expand.grid(cutoffQ_num, maxGap_int)

  for(i in 1:nrow(designPts_mat)){

    ###  Generate Data Set  ###
    delta <- designPts_mat[i, 1]
    seed  <- designPts_mat[i, 2]

    treatment_ls <- SimulateData(beta_mat = beta_mat,
                                 AclustCPG_df = Aclusters_df,
                                 delta_num = delta,
                                 seed_int = seed)
    betas_df <- treatment_ls$simBetaVals_df


    ###  Data Wrangling  ###
    mergedBetas_df <- merge(betas_df, CPGs_df,
                            by.x = "row.names",
                            by.y = "ILMNID")

    cpgInfo_df <- subset(mergedBetas_df, select = c("chr", "MAPINFO"))
    cpgInfo_df$chr <- substr(cpgInfo_df$chr, 4, 6)

    betaSorted_df <- mergedBetas_df
    row.names(betaSorted_df) <- betaSorted_df$Row.names
    betaSorted_df$Row.names <-
      betaSorted_df$chr <-
      betaSorted_df$MAPINFO <-
      NULL
    betaSorted_mat <- as.matrix(betaSorted_df)


    ###  Inner Parameter Grid Search  ###
    for(j in 1:nrow(paramsGrid_mat)){

      ###  Calculate Method Output  ###
      cutoffQ <- paramsGrid_mat[j, 1]
      maxGap  <- paramsGrid_mat[j, 2]

      ###  Parallel Setup  ###
      if(!parallel){
        numCores <- 1
      }
      # Clean memory
      rm(treatment_ls, betas_df, mergedBetas_df, betaSorted_df)
      # Make and Register Cluster
      clust <- makeCluster(numCores)
      registerDoParallel(clust)

      suppressMessages(
        res_ls <- RunBumphunter(betaVals_mat = betaSorted_mat,
                                chromos_char = cpgInfo_df$chr,
                                chromPosit_num = cpgInfo_df$MAPINFO,
                                cpgLocation_df = CPGs_df,
                                pickCutoffQ_num = cutoffQ,
                                maxGap_int = maxGap,
                                numCores = numCores)
      )

      stopCluster(clust)

      ###  Define NULL Data  ###
      if(is.null(res_ls[[1]])){
        res_ls[[1]] <- data.frame(
          dmr.chr     = NA_character_,
          dmr.start   = NA_integer_,
          dmr.end     = NA_integer_,
          chr         = NA_character_,
          start       = NA_integer_,
          end         = NA_integer_,
          value       = NA_real_,
          area        = NA_real_,
          cluster     = NA_integer_,
          indexStart  = NA_integer_,
          indexEnd    = NA_integer_,
          L           = NA_integer_,
          clusterL    = NA_integer_,
          p.value     = NA_real_,
          fwer        = NA_real_,
          p.valueArea = NA_real_,
          fwerArea    = NA_real_,
          dmr.pval    = NA_real_,
          dmr.n.cpgs  = NA_integer_
        )
      }

      ###  Save Results  ###
      file_char <- paste0(
        resultsDir, "BumphunterResults_delta", delta, "_seed", seed,
        "_pickQ", cutoffQ, "_maxGap", maxGap, ".RDS"
      )

      if(verbose){
        message("Saving results to file ", file_char, "\n")
      }

      saveRDS(res_ls, file = file_char)

    } # END for(j)

  } # END for(i)

}
