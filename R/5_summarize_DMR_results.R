#' Summarize Method Accuracy of Returned Results
#'
#' @description Given the output of the \code{\link{CleanResults}} (cleaned and
#'    joined predicted clusters to true clusters), return a data frame of the
#'    precision, power, elapsed time, and other summary statistics.
#'
#' @param cleanDMR_df A data frame of joined predicted and true clusters, as
#'    returned by the \code{\link{CleanResults}} function.
#' @param time_num The computing time for the given method and parameter
#'    configuration, as stored in the second entry of the output list from the
#'    \code{\link{WriteDMRcateResults}}, \code{\link{WriteBumphunterResults}},
#'    or \code{\link{WriteProbeLassoResults}} functions.
#'
#' @return A data frame of model fit statistics for the given method to the
#'    data generated according to the given parameter configuration:
#'    \itemize{
#'      \item{\code{time}}{},
#'      \item{\code{FN}}{},
#'      \item{\code{FP}}{},
#'      \item{\code{TN}}{},
#'      \item{\code{TP}}{},
#'      \item{\code{power}}{},
#'      \item{\code{nPower}}{},
#'      \item{\code{FPprecis}}{},
#'      \item{\code{TPprecis}}{},
#'      \item{\code{precision}}{},
#'      \item{\code{nPrecis}}{},
#'      \item{\code{nCPG_q1}}{},
#'      \item{\code{nCPG_med}}{},
#'      \item{\code{nCPG_q3}}{}
#'    }
#'
#' @importFrom PRROC pr.curve
#' @importFrom stats aggregate
#' @importFrom stats as.formula
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called internally by the ProcessDMRcateResults(), ProcessCombpResults(),
#'   #   ProcessProbeLassoResults(), and ProcessBumphunterResults() functions.
SummarizeResults <- function(cleanDMR_df, time_num){
  # browser()

  ###  Table Power Results  ###
  # Frequency count of each status - based on unique aclusters, for power
  #   calculation
  statusCls_df <- unique(cleanDMR_df[, c("Clusternumber", "status")])

  powerSummary_tbl <- table(
    factor(statusCls_df$status,
           levels = c("FN", "FP", "TN", "TP"))
  )
  powerSummary_df <- data.frame(
    matrix(as.numeric(powerSummary_tbl), ncol = 4)
  )
  colnames(powerSummary_df) <- names(powerSummary_tbl)


  ###  Add Power Calculation Summary ###
  powerSummary_df$power <-
    powerSummary_df$TP / (powerSummary_df$TP + powerSummary_df$FN)
  powerSummary_df$time <- time_num
  powerSummary_df$nPower <- powerSummary_df$TP + powerSummary_df$FN
  powerSummary_df <-
    powerSummary_df[, c("time", "FN", "FP", "TN", "TP", "power", "nPower")]


  ###  Area under Precision-Recall Curve  ###
  if(!is.null(cleanDMR_df$dmr.pval)){

    x_df <- cleanDMR_df[, c("aclust.order", "dmr.pval", "actual")]
    x_df$status <- ifelse(x_df$actual == "positive", 1, 0)
    x_df$dmr.pval[is.na(x_df$dmr.pval)] <- 1

    # take min pvalue of dmrs
    agg_fmla <- as.formula("dmr.pval ~ aclust.order + actual + status")
    x_df <- aggregate(agg_fmla, data = x_df, FUN = min)

    # PR curve
    prCurve_ls <- pr.curve(scores.class0 = 1 - x_df$dmr.pval,
                           weights.class0 = x_df$status)
    powerSummary_df$AuPR <- prCurve_ls$auc.integral
    rm(x_df, agg_fmla, prCurve_ls)

  } else {
    powerSummary_df$AuPR <- NA_real_
  }


  ###  Table Precision Results  ###
  # Frequency count of each status - based on unique DMRs, for precision
  #   calculation
  if(!is.null(cleanDMR_df$dmr.order)){

    statusDMR_df <- unique(cleanDMR_df[, c("dmr.order", "status")])

    precisSummary_tbl <- table(
      factor(statusDMR_df$status,
             levels = c("FN", "FP", "TN", "TP"))
    )
    precisSummary_df <- data.frame(
      matrix(as.numeric(precisSummary_tbl), ncol = 4)
    )
    colnames(precisSummary_df) <- paste0(names(precisSummary_tbl), "precis")
    precisSummary_df$FNprecis <- precisSummary_df$TNprecis <- NULL


    ###  Add Precision Calculation Summary ###
    precisSummary_df$precision <-
      precisSummary_df$TP / (precisSummary_df$TP + precisSummary_df$FP)
    precisSummary_df$nPrecis <- precisSummary_df$TP + precisSummary_df$FP


    powerSummary_df <- cbind(powerSummary_df, precisSummary_df)

  } else {

    powerSummary_df$FPprecis  <- NA_integer_
    powerSummary_df$TPprecis  <- NA_integer_
    powerSummary_df$precision <- NA_real_
    powerSummary_df$nPrecis   <- NA_integer_

  }


  ###  Precision Statistics  ###
  # Matthews Correlation Coefficient
  CalcMCC <- function(tp, tn, fp, fn){
    (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  }
  powerSummary_df$mcc <- CalcMCC(tp = powerSummary_df$TPprecis,
                                 tn = powerSummary_df$TN,
                                 fp = powerSummary_df$FPprecis,
                                 fn = powerSummary_df$FN)

  # F1-score
  CalcF1 <- function(tp, fp, fn){ (2 * tp) / (2 * tp + fp + fn) }
  powerSummary_df$F1  <- CalcF1(tp = powerSummary_df$TPprecis,
                                fp = powerSummary_df$FPprecis,
                                fn = powerSummary_df$FN)


  ###  Summary of numCPGs  ###
  # I've tried to fit the Poisson and Gamma distributions to the number of CPGs,
  #   and a Gamma to the log of the number of CPGs. Nothing fits well. The
  #   counts are minima-inflated. We will report the three quartiles.
  nCPGsSummary_num <- summary(cleanDMR_df$dmr.n.cpgs)
  nCPGs_df <- data.frame(
    matrix(nCPGsSummary_num[c(2, 3, 5)], nrow = 1),
    stringsAsFactors = FALSE
  )
  colnames(nCPGs_df) <- paste0("nCPG_", c("q1", "med", "q3"))

  powerSummary_df <- cbind(powerSummary_df, nCPGs_df)


  ###  Return  ###
  powerSummary_df

}
