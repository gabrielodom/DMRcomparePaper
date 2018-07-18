#' Merge CPGs with DMRs
#'
#' @description An internal wrangling function for \code{RunDMRcate}
#'
#' @param DMRs_df A data frame with variables \code{dmr.start}, \code{dmr.end},
#'    \code{dmr.chr} (e.g. \code{chr1}), and \code{dmr.pval}. This data frame
#'    will be returned by the \code{extractRanges()} function from the
#'    \code{DMRcate} package.
#' @param CPGs_df A data frame matching chromosomes to CPG names and locations.
#'    This data frame contains the variables \code{ILMNID}, \code{MAPINFO}, and
#'    \code{chr} (e.g. \code{chr1}) and is given in the \code{cpgLocation_df}
#'    data set.
#' @param alpha Significance level to select regions (with \code{dmr.pval} less
#'    than the specified value). Defaults to 0.05.
#'
#' @return A data frame containing DMR info, CPG info, and the number of CPGs
#'
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges findOverlaps
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called internally by the StandardizeOutput() function.
MergeDMRsWithCPGs <- function(DMRs_df, CPGs_df, alpha = 0.05){

  ###  1. Make ranges of DMR info  ###
  DMRs_df <- DMRs_df[DMRs_df$dmr.pval < alpha , ]
  sig.ranges <- IRanges(DMRs_df$dmr.start, DMRs_df$dmr.end)
  dmr.ranges <- GRanges(seqnames = DMRs_df$dmr.chr, ranges = sig.ranges)


  ###  2. Make ranges of CPG info  ###
  temp.cpg.ranges <- IRanges(CPGs_df$MAPINFO, CPGs_df$MAPINFO)
  cpg.ranges <- GRanges(seqnames = CPGs_df$chr, ranges = temp.cpg.ranges)


  ###  3. Find overlaps  ###
  overlaps_df <- as.data.frame(
    findOverlaps(query = dmr.ranges, subject = cpg.ranges, type = "any")
  )

  overlaps_df$dmr.order <- as.numeric(overlaps_df$queryHits)
  overlaps_df$cpg.order <- as.numeric(overlaps_df$subjectHits)


  ###  4. Merge with DMR info  ###
  DMRs_df$dmr.row <- 1:nrow(DMRs_df)
  DMRsInfo_df <- merge(x = overlaps_df, y = DMRs_df,
                       by.x = "dmr.order", by.y = "dmr.row")


  ###  5. Merge with CPGs info  ###
  CPGs_df$row <- 1:nrow(CPGs_df)
  overlapInfo_df <- merge(DMRsInfo_df, CPGs_df,
                          by.x = "cpg.order", by.y = "row")
  overlapInfo_df <- overlapInfo_df[order(overlapInfo_df$dmr.order), ]


  ###  6. Add number of CPGs and Return  ###
  numCPGs <- as.data.frame(table(overlapInfo_df$dmr.order))

  merge(overlapInfo_df, numCPGs,
        by.x = "dmr.order", by.y = "Var1")

}
