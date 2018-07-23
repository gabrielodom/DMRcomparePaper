#' Annotated CPG Data Set
#'
#' @description Trimmed beta values for 14 subjects over a 450k methylation
#'    array.
#'
#' @format A matrix containing 356603 CPGs measured on 14 subjects. These
#'    values have been trimmed to only range from 0.05 to 0.95. The column
#'    names indicate from which phenotypic group the subjects were drawn.
#'
#' @keywords internal
#'
#' @source Calculated via the \code{1_Aclust_data_import.R} script in the
#'    \code{old_scripts} sub-directory of the \code{inst} directory.
"betaVals_mat"
