#' Annotated CPG Data Set
#'
#' @description A beta value matrix for selected methylation samples from a 450k
#'    methylation array. These beta values are trimmed to range from 0.05 to
#'    0.95.
#'
#' @format A matrix containing 356603 CPGs measured on 14 subjects. The rownames
#'    are the CPG IDs, and the column names indicate the sample IDs of the
#'    subjects.
#'
#' @source Calculated via the \code{1_Aclust_data_import.R} script in the
#'    \code{old_scripts} sub-directory of the \code{inst} directory.
"betaVals_mat"
