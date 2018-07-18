#' CPG List-to-Table Conversion
#'
#' @description Convert a list of CPGs to a table
#'
#' @param cpgs_ls a list where each item has CPG IDs in a cluster
#' @param methylval_df a dataframe of beta values
#'
#' @return a data table with
#'    \itemize{
#'      \item{\code{rows}}{cluster 1,1,1, 2,2,2}
#'      \item{\code{columns}}{sample IDs}
#'      \item{\code{value}}{mvalue or beta values}
#'    }
#'
#' @importFrom data.table rbindlist
#'
#' @export
#'
#' @keywords internal
#'
#' @examples
#'   # Called by the script in inst/old_scripts/1_Aclust_data_import.R
ConvertCPGList <- function(cpgs_ls, methylval_df){

  cluster_ls <- lapply(cpgs_ls, function(item){
    as.data.frame(methylval_df[item, ])
  })
  clusterRowLabel <- lapply(cpgs_ls, function(item){
    as.data.frame(rownames(methylval_df[item, ]))
  })

  cluster_tab <- rbindlist(cluster_ls,
                           idcol = "cluster",
                           use.names = TRUE)
  clusterLabel_tab <- rbindlist(clusterRowLabel,
                                idcol = "cluster",
                                use.names = TRUE)
  colnames(clusterLabel_tab)[2] <- "probeID"

  cbind(clusterLabel_tab[, 2], cluster_tab)
}
