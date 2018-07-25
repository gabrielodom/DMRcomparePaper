#' CPG Locations
#'
#' @description An annotation table that indicates locations of CpGs. This data
#'    frame has CPG IDs as the rows with matching chromosome and location info
#'    in the columns.
#'
#' @format A data frame containing 485512 CPGs (rows) and three columns. The
#'    columns are:
#'    \itemize{
#'      \item{\code{ILMNID} : }{the CPG ID, as a factor}
#'      \item{\code{chr} : }{the chromosome label, as a character}
#'      \item{\code{MAPINFO} : }{the chromosome location, as an integer}
#'    }
#'
#' @source Compiled via the \code{1_Aclust_data_import.R} script in the
#'    \code{old_scripts} sub-directory of the \code{inst} directory.
"cpgLocation_df"
