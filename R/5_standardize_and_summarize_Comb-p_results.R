#' Extract and Process Comb-p Results Files
#'
#' @description Given a directory of saved Comb-p results, as \code{.RDS} files,
#'    import, standardize, and summarize these data files.
#'
#' @param resultsDir The name of the directory where the Comb-p method results
#'    are stored.
#' @param beta_mat A matrix of beta values across genome on the array. This is
#'    given in the \code{betaVals_mat} data set.
#' @param AclustCPG_df A data frame of \code{Aclust} results. This is given in
#'    the \code{startEndCPG_df} data set.
#' @param cpgLocation_df A data frame matching chromosomes to CPG names and
#'    locations. This is given in the \code{cpgLocation_df} data set.
#' @param dmr.sig.threshold Significance level to select regions (with
#'    \code{dmr.pval} less than the specified value) passed to the internal
#'    \code{\link{MergeDMRsWithCPGs}} function.
#' @param min.cpgs The minimum number of CPGs necessary to consider a result
#'    significant. Defaults to 5.
#' @param verbose Should the function print progress messages? Defaults to
#'    \code{TRUE}.
#'
#' @return A data frame of model fit statistics for the Comb-p method under
#'    each of the given parameter combinations to the data generated for each
#'    design configuration
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   data("betaVals_mat")
#'   data("startEndCPG_df")
#'   data("cpgLocation_df")
#'
#'   combpRes_df <- ProcessCombpResults(
#'     resultsDir = "DMRcate_results/",
#'     beta_mat = betaVals_mat,
#'     AclustCPG_df = startEndCPG_df,
#'     cpgLocation_df = cpgLocation_df
#'   )
#' }
ProcessCombpResults <- function(resultsDir,
                                beta_mat,
                                AclustCPG_df,
                                cpgLocation_df,
                                dmr.sig.threshold = 0.05,
                                min.cpgs = 5,
                                verbose = TRUE){
  # browser()

  ###  Vector of Appropriate File Names  ###
  files_char <- list.files(path = resultsDir, pattern = "^Combp.*RDS$")
  files_char <- gsub(pattern = ".RDS", replacement = "", files_char)
  splitFileNames_ls <- strsplit(files_char, split = "_")


  ###  Initialize Design Points  ###
  ExtractParamLevels <- function(splitNames_ls, param){

    posIdx <- grep(param, splitNames_ls[[1]])
    point_char <- sapply(splitNames_ls, `[`, posIdx)
    point_num <- as.numeric(
      gsub(pattern = param, replacement = "", point_char)
    )

    sort(unique(point_num))

  }

  deltas_num <- ExtractParamLevels(splitFileNames_ls, param = "delta")
  seeds_int <- ExtractParamLevels(splitFileNames_ls, param = "seed")
  combSeeds_num <- ExtractParamLevels(splitFileNames_ls, param = "combSeed")
  combDists_int <- ExtractParamLevels(splitFileNames_ls, param = "combDist")

  # This ensures that the delta values stay together, not the seed values.
  designPts_mat <- expand.grid(seeds_int, deltas_num)
  paramsGrid_mat <- expand.grid(combSeeds_num, combDists_int)


  ###  Simulate Gold Standard Outer Loop  ###
  out_ls <- vector(mode = "list", length = nrow(designPts_mat))
  for(i in 1:nrow(designPts_mat)){
    # browser()

    ###  Generate Data Set  ###
    seed  <- designPts_mat[i, 1]
    delta <- designPts_mat[i, 2]

    treatment_ls <- SimulateData(beta_mat = beta_mat,
                                 AclustCPG_df = AclustCPG_df,
                                 delta_num = delta,
                                 seed_int = seed)
    trueClusters_df <- treatment_ls$simAclusters_df


    ###  Inner Results Comparison  ###
    innerOut_ls <- vector(mode = "list", length = nrow(paramsGrid_mat))
    for(j in 1:nrow(paramsGrid_mat)){
      # browser()

      ###  Calculate Method Output  ###
      cSeed <- paramsGrid_mat[j, 1]
      cDist <- paramsGrid_mat[j, 2]


      ###  Load Results  ###
      resFileName_char <- paste0(resultsDir, "CombpResults",
                                 "_delta", delta, "_seed", seed,
                                 "_combSeed", cSeed, "_combDist", cDist,
                                 ".RDS")
      res_ls <- readRDS(resFileName_char)


      ###  Extract and Standardize Output  ###
      # Remove any rows with all NAs
      ranges_df <- res_ls[[1]]
      ranges_df <- ranges_df[rowSums(is.na(ranges_df)) != ncol(ranges_df), ]
      # Standardize: if all of the rows of the results df are removed, then
      #   skip straight to the "clean results" step.
      if(nrow(ranges_df) > 0){

        res_ls[[1]] <- StandardizeOutput(
          methodOut_df = ranges_df,
          method = "Comb_p",
          cpgLocation_df = cpgLocation_df,
          dmr.sig.threshold = dmr.sig.threshold,
          min.cpgs = min.cpgs
        )

      }


      ###  Clean and Summarize Results  ###
      all_df <- CleanResults(dmrResults_ls = res_ls,
                             Aclusters_df = trueClusters_df)
      innerOut_df <- SummarizeResults(cleanDMR_df = all_df,
                                      time_num = res_ls[[2]])
      innerMeta_df <- data.frame(method = "Comb-p",
                                 delta = delta,
                                 seed = seed,
                                 combSeed = cSeed,
                                 combDist = cDist,
                                 stringsAsFactors = FALSE)
      innerOut_ls[[j]] <- cbind(innerMeta_df, innerOut_df)

    } # END for(j)


    ###  Bind Inner Results List  ###
    out_ls[[i]] <- do.call(rbind, innerOut_ls)

    if(verbose){
      message("Completed summary for delta = ", delta,
              " and seed = ", seed, ".")
    }


  } # END for(i)

  ###  Bind Outer Results List and Return  ###
  out_df <- do.call(rbind, out_ls)
  # This would return the results to the console, but we should probably save
  #   the results data frame instead?
  out_df

}
