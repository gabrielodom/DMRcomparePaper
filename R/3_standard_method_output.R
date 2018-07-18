#' Standardize DMR Method Comparison Output
#'
#' @description Given the output from either \code{\link[DMRcate]{dmrcate}},
#'    \code{\link[ChAMP]{champ.DMR}} (for the ProbeLasso method), or
#'    \code{\link[bumphunter]{bumphunter}}, standardize the $p$-value, start
#'    and end of each DMR, and chromosome columns.
#'
#' @param methodOut_df An object returned by \code{\link[DMRcate]{dmrcate}},
#'    \code{\link[ChAMP]{champ.DMR}} (with \code{method = "ProbeLasso"}), or
#'    \code{\link[bumphunter]{bumphunter}}.
#' @param method A character string matching \code{"DMRcate"},
#'    \code{"ProbeLasso"}, or \code{"Bumphunter"}. Partial matching is handled
#'    via the \code{\link[base]{match.arg}} function.
#' @param cpgLocation_df A data frame matching chromosomes to CPG names and
#'    locations. This is given in the \code{cpgLocation_df} data set.
#' @param dmr.sig.threshold Significance level to select regions (with
#'    \code{dmr.pval} less than the specified value) passed to the internal
#'    \code{\link{MergeDMRsWithCPGs}} function.
#' @param min.cpgs The minimum number of CPGs before we consider a result
#'    significant. Defaults to 5.
#'
#' @return A data frame of ...
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called internally by the RunDMRcate(), RunProbeLasso(), RunBumphunter(),
#'   #   ProcessCombpResults(), BuildOverlaps(), and BuildPRcurve() functions.
StandardizeOutput <- function(methodOut_df,
                              method = c("DMRcate",
                                         "ProbeLasso",
                                         "Bumphunter",
                                         "Comb_p"),
                              cpgLocation_df,
                              dmr.sig.threshold = 0.05,
                              min.cpgs = 5){

  method <- match.arg(method)

  ###  Extract Output  ###
  switch(method,

    DMRcate = {

      methodOut_df$dmr.pval  <- methodOut_df$Stouffer
      methodOut_df$dmr.chr   <- methodOut_df$seqnames
      methodOut_df$dmr.start <- methodOut_df$start
      methodOut_df$dmr.end   <- methodOut_df$end

    },
    ProbeLasso = {

      methodOut_df$dmr.chr   <- methodOut_df$seqnames
      methodOut_df$dmr.start <- methodOut_df$start
      methodOut_df$dmr.end   <- methodOut_df$end
      methodOut_df$dmr.pval  <- methodOut_df$dmrP

    },
    Bumphunter = {

      methodOut_df$dmr.pval  <- methodOut_df$p.valueArea
      methodOut_df$dmr.chr   <- paste0("chr", as.character(methodOut_df$chr))
      methodOut_df$dmr.start <- methodOut_df$start
      methodOut_df$dmr.end   <- methodOut_df$end

    },
    Comb_p = {

      methodOut_df$dmr.pval  <- methodOut_df$z_p
      methodOut_df$dmr.chr   <- methodOut_df$chrom
      methodOut_df$dmr.start <- methodOut_df$start
      methodOut_df$dmr.end   <- methodOut_df$end

    }

  )


  ###  Transform  ###
  temp <- MergeDMRsWithCPGs(DMRs_df = methodOut_df,
                            CPGs_df = cpgLocation_df,
                            alpha = dmr.sig.threshold)
  temp.ncpgs <- unique(
    subset(temp, select = c("dmr.chr", "dmr.start", "dmr.end", "Freq"))
  )

  clean_df <-
    merge(methodOut_df, temp.ncpgs, by = c("dmr.chr", "dmr.start", "dmr.end"))
  clean_df$dmr.n.cpgs <- clean_df$Freq
  clean_df$Freq <- NULL
  clean_df <- clean_df[clean_df$dmr.n.cpgs >= min.cpgs, ]

  if(nrow(clean_df) == 0){
    clean_df <- NULL
  }


  ###  Return  ###
  clean_df

}
