#' Clean and Join DMR Method Predictions to True Values
#'
#' @description Given a results file, as written by one of the
#'    \code{Write*Results} functions, and the data frame of true cluster values,
#'    as returned by the \code{SimulateData}, join the predicted clusters to the
#'    true clusters.
#'
#' @param dmrResults_ls A list of results, written by one of the functions
#'    \code{\link{WriteDMRcateResults}}, \code{\link{WriteBumphunterResults}},
#'    or \code{\link{WriteProbeLassoResults}}.
#' @param Aclusters_df A data frame of true cluster values, returned in the
#'    \code{simAclusters} entry of the output of the \code{\link{SimulateData}}
#'    function.
#'
#' @return A data frame with the predicted clusters joined to their true values.
#'    This data frame is to be immediately passed to the
#'    \code{\link{SummarizeResults}} function for analysis.
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called internally by the ProcessDMRcateResults(), ProcessCombpResults(),
#'   #   ProcessProbeLassoResults(), ProcessBumphunterResults(), and
#'   #   BuildPRcurve() functions.
CleanResults <- function(dmrResults_ls, Aclusters_df) {

  ranges_df <- dmrResults_ls[[1]]
  # Remove rows with all NAs
  ranges_df <- ranges_df[rowSums(is.na(ranges_df)) != ncol(ranges_df), ]
  elapsedtime <- dmrResults_ls[[2]]

  clusters_df <- unique(
    Aclusters_df[, c("Clusternumber",
                     "chromosome",
                     "start_position",
                     "end_position",
                     "actual")]
  )

  if(nrow(ranges_df) > 0){

    ###  Create GRanges  ###
    # query = significant DMRs; need to limit to min.cpgs > 4 and pval < 0.05
    signifRanges_df <-
      ranges_df[ranges_df$dmr.n.cpgs > 4 & ranges_df$dmr.pval < 0.05, ]
    query_GR <- GRanges(seqnames = signifRanges_df$dmr.chr,
                        ranges = IRanges(signifRanges_df$dmr.start,
                                         signifRanges_df$dmr.end))

    # subject = Aclusters
    subject_GR <- GRanges(seqnames = clusters_df$chromosome,
                          ranges = IRanges(clusters_df$start_position,
                                           clusters_df$end_position))

    # subject-query overlap
    overlap_df <- as.data.frame(
      findOverlaps(query_GR, subject_GR, type = "any", select = "all")
    )
    overlap_df$dmr.order <- as.numeric(overlap_df$queryHits)
    overlap_df$aclust.order <- as.numeric(overlap_df$subjectHits)


    ###  Merge Results  ###
    # merge with dmr info
    signifRanges_df$dmr.row <- 1:nrow(signifRanges_df)
    signifRanges_df$predicted <- "positive"
    overlapDMRs_df <- merge(x = overlap_df, y = signifRanges_df,
                            by.x = "dmr.order", by.y = "dmr.row", all = TRUE)

    # merge with aclust info
    clusters_df$aclust.row <- 1:nrow(clusters_df)

    # merge all
    all_df <- merge(x = overlapDMRs_df, y = clusters_df,
                    by.x = "aclust.order", by.y = "aclust.row", all = TRUE)

  } else {

    all_df <- clusters_df
    all_df$predicted <- "negative"

  }

  ###  Add Status Column  ###
  all_df$predicted[is.na(all_df$predicted)] <- "negative"
  all_df$actual[is.na(all_df$actual)] <- "negative"

  # True Positive
  all_df$status[
    all_df$actual == "positive" & all_df$predicted == "positive"
    ] <- "TP"

  # False Negative
  all_df$status[
    all_df$actual == "positive" & all_df$predicted == "negative"
    ] <- "FN"

  # False Positive
  all_df$status[
    all_df$actual == "negative" & all_df$predicted == "positive"
    ] <- "FP"

  # True Negative
  all_df$status[
    all_df$actual == "negative" & all_df$predicted == "negative"
    ] <- "TN"

  ###  Return  ###
  all_df

}

