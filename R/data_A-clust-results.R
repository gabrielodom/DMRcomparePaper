#' Annotated CPG Data Set
#'
#' @description Beta values for 14 subjects over a 450k methylation array.
#'
#' @format A data frame containing 20361 CPG locations measured on 14 subjects.
#'    The first eight columns are the metadata for the CPG, including starting
#'    and ending positions on the chromosome. The column names for the subjects
#'    indicate from which phenotypic group the subjects were drawn. The row
#'    names are the CPGs.
#'
#' @keywords internal
#'
#' @source Calculated via the \code{1_Aclust_data_import.R} script in the
#'    \code{old_scripts} sub-directory of the \code{inst} directory.
"startEndCPG_df"
