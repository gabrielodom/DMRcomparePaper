#' Annotated CPG Data Set
#'
#' @description A matrix of beta values for clusters of CpGs over a 450k
#'    methylation array.
#'
#' @format A data frame containing 20361 CPG locations measured on 14 subjects.
#'    The rows are the CPG IDs. The first eight columns are the metadata for
#'    the CPGs, including: \code{Clusternumber}, the CPG ID (\code{cpg}),
#'    chromosome (\code{CHR}), chromosome location (\code{MAPINFO}), chromosome
#'    start position (\code{start_position}), and chromosome end position
#'    (\code{end_position}). The remaining 14 columns are the beta values for
#'    the subjects. The column names for the subjects indicate from which
#'    phenotypic group the subjects were drawn; for example, the
#'    \code{9744-Tumor} column indicates that this subject was from the case
#'    group.
#'
#' @source Calculated via the \code{1_Aclust_data_import.R} script in the
#'    \code{old_scripts} sub-directory of the \code{inst} directory.
"startEndCPG_df"
