#' Build a List of Precision-Recall Curve Objects
#'
#' @description Given a directory of best-performing results files from one of
#'    the simulation functions (\code{\link{WriteDMRcateResults}},
#'    \code{\link{WriteProbeLassoResults}},
#'    \code{\link{WriteBumphunterResults}}, or results from the \code{Comb-p}
#'    method in \code{Python}), import the raw data files and construct PR-curve
#'    objects via the \code{\link[PRROC]{pr.curve}} function.
#'
#' @param bestResultsDir The name of the directory where the method results from
#'    the best-performing parameter settings are stored. For the full design we
#'    have included (\code{delta = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)} and
#'    \code{seed = c(100, 210, 330, 450, 680)}), this directory should contain
#'    35 \code{.RDS} files per method.
#' @param delta A treatment size corresponding to one of the simulations with
#'    completed results files in the \code{bestResultsDir} directory.
#' @param seed A seed value corresponding to one of the simulations with
#'    completed results files in the \code{bestResultsDir} directory.
#' @param beta_mat A matrix of beta values across genome on the array. The
#'    default value is given in the \code{betaVals_mat} data set.
#' @param AclustCPG_df A data frame of \code{Aclust} results. The default value
#'    is given in the \code{startEndCPG_df} data set.
#' @param CPGs_df A data frame matching chromosomes to CPG names and
#'    locations. The default value is given in the \code{cpgLocation_df} data
#'    set, passed to the \code{\link{StandardizeOutput}} function. This data
#'    set is only necessary if the results directory contains Comb-p results
#'    with the specified \code{delta} and \code{seed} values.
#' @param min.cpgs The minimum number of CPGs necessary to consider a result
#'    significant. Defaults to 5. This argument is only required if the results
#'    directory contains Comb-p results with the specified \code{delta} and
#'    \code{seed} values.
#'
#' @return A list of PR-curve objects, to be plotted via the
#'    \code{\link{PlotPRCurve}} function.
#'
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#' @importFrom PRROC pr.curve
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    BuildPRcurve(
#'      bestResultsDir = "best_cases_results/",
#'      delta = 0.4,
#'      seed = 100
#'    )
#' }
BuildPRcurve <- function(bestResultsDir,
                         delta = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4),
                         seed = c(100, 210, 330, 450, 680),
                         beta_mat = betaVals_mat,
                         AclustCPG_df = startEndCPG_df,
                         CPGs_df = cpgLocation_df,
                         min.cpgs = 5){

  ###  Gold Standard  ###
  treatment_ls <- SimulateData(beta_mat = beta_mat,
                               AclustCPG_df = AclustCPG_df,
                               delta_num = delta,
                               seed_int = seed)
  trueClusters_df <- treatment_ls$simAclusters_df

  ### List Results Files  ###
  fileNames_char <- list.files(bestResultsDir)
  targetNames_char <- paste0("delta", delta, "_seed", seed)
  correctFiles_idx <- grep(targetNames_char, fileNames_char)
  correctNames_char <- fileNames_char[correctFiles_idx]
  names(correctNames_char) <- sapply(
    strsplit(correctNames_char, "_"),
    function(x){
      gsub(pattern = "Results", replacement = "", x[1])
    }
  )

  ### Load and Clean Appropriate Files  ###
  allRes_ls <- lapply(correctNames_char, function(x){

    res_ls <- readRDS(paste0(bestResultsDir, x))

    if(grepl("Combp", x)){

      results_df <- StandardizeOutput(
        methodOut_df = res_ls[[1]],
        method = "Comb_p",
        cpgLocation_df = CPGs_df
      )

      # The raw Comb-p results were filtered to > 1, not > 4.
      results_df <- results_df[results_df$dmr.n.cpgs >= min.cpgs, ]
      res_ls[[1]] <- results_df

    }

    CleanResults(dmrResults_ls = res_ls,
                 Aclusters_df = trueClusters_df)

  })

  ###  Build PR Curves  ###
  allPRs_ls <- lapply(allRes_ls, function(cleanDMR_df){

    if(!is.null(cleanDMR_df$dmr.n.cpgs)){
      # This is because ProbeLasso randomly sucks for small delta.

      x_df <- cleanDMR_df[, c("aclust.order", "dmr.pval", "actual")]
      x_df$status <- ifelse(x_df$actual == "positive", 1, 0)
      x_df$dmr.pval[is.na(x_df$dmr.pval)] <- 1

      # take min pvalue of dmrs
      agg_fmla <- as.formula("dmr.pval ~ aclust.order + actual + status")
      x_df <- aggregate(agg_fmla, data = x_df, FUN = min)

      # PR curve
      pr.curve(scores.class0 = 1 - x_df$dmr.pval,
               weights.class0 = x_df$status,
               curve = TRUE)

    } else {
      NULL
    }

  })

  ###  Return  ###
  attr(allPRs_ls, "delta") <- delta
  attr(allPRs_ls, "repl")  <- which(c(100, 210, 330, 450, 680) %in% seed)
  allPRs_ls

}
