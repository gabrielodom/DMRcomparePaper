#' Extract and Process Comb-p Results Files
#'
#' @description Given a directory of saved Comb-p results, as \code{.RDS} files,
#'    import, standardize, and summarize these data files.
#'
#' @param resultsDir The name of the directory where the Comb-p method results
#'    are stored.
#'    
#' @param beta_mat A beta value matrix for methylation samples from a
#'    450k methylation array with CpG IDs as the row names and sample IDs as
#'    the column names. An example is given in the \code{betaVals_mat} data set.
#' 
#' @param AclustCPG_df A data frame of beta values and CpG information for
#'    clusters of CpGs over a 450k methylation array. The rows correspond to the
#'    CPGs. The columns have information on the cluster number, chromosome,
#'    cluster start and end locations, and the beta values for each subject
#'    grouped by some clinical indicator (e.g. case v. control). An example is
#'    given in the \code{startEndCPG_df} data set. This data set also has
#'    information on true status of the clusters, via variable \code{status}, 
#'    with values "positive" or "negative", indicating whether treatment effect
#'    was added to the cluster.
#' 
#' @param cpgLocation_df An annotation table that indicates locations of CpGs.
#'    This data frame has CpG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set.
#'    
#' @param dmr.sig.threshold Regions with DMR p-value less than
#'    \code{dmr.sig.threshold} are selected for the output. 
#' 
#' @param min.cpgs Minimum number of CpGs. Regions with at least \code{min.cpgs}
#'    are selected for the output. Defaults to 5.
#'    
#' @param verbose Should the function print progress messages? Defaults to
#'    \code{TRUE}.
#'
#' @return A data frame of model performance measures for the Comb-p method
#'    under each of the given parameter combinations applied to the data
#'    generated with different treatment effects 
#'    
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
