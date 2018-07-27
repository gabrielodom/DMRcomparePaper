#' Simulate Differences in Methylation Data
#'
#' @description This function generates a methylation dataset where treatment effects (some constant values) are added 
#' to to beta values in one observation class for some randomly selected co-methylated clusters 
#'    
#'
#' @param beta_mat A beta value matrix for methylation samples from a
#'    450k methylation array with CPG IDs as the row names and sample IDs as
#'    the column names. An example is given in the \code{betaVals_mat} data set.
#' @param AclustCPG_df A data frame of beta values and CpG information for
#'    clusters of CpGs over a 450k methylation array. The rows correspond to the
#'    CPGs. The columns have information on the cluster number, chromosome,
#'    cluster start and end locations, and the beta values for each subject
#'    grouped by some clinical indicator (e.g. case v. control). An example is
#'    given in the \code{startEndCPG_df} data set.
#' @param delta_num The treatment size: a non-negative real number to add to
#'    the beta values within randomly-selected clusters for a single class of
#'    subjects. This artifically creates differentially-methylated regions
#'    (DMRs).
#' @param seed_int The seed value passed to the \code{\link[base]{Random}}
#'    function to enable reproducible results
#' @param betaCols_idx The column numbers of the \code{AclustCPG_df} data frame
#'    in which beta values for each subject are stored. This function assumes
#'    that the subject columns are grouped by their class.
#' @param numEx_int The number of samples in the first group. Once again, this
#'    function assumes that these samples are contiguous columns of the
#'    \code{AclustCPG_df} data frame.
#' @param numClusters_int The total number of clusters to randomly select to
#'    be inflated by the treatment amount, \code{delta_num}
#'
#' @return A list with two elements:
#'    \itemize{
#'      \item{\code{simBetaVals_df}}{A data frame of beta values after
#'        treatment effects were added, used for input for different DMR-finding
#'        methods. Note this is whole-genome data.}
#'      \item{\code{simAclusters_df}}{A data frame of the methylation values
#'        only for \code{Aclust} and annotation for whether treatment effects
#'        were added. Note this has only CPGs mapped to all the clusters found
#'        by the \code{Aclust} method.}
#'    }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    data("startEndCPG_df")
#'    data("betaVals_mat")
#'
#'    SimulateData(beta_mat = betaVals_mat,
#'                 AclustCPG_df = startEndCPG_df,
#'                 delta_num = 0.4,
#'                 seed_int = 12345)
#' }
SimulateData <- function(beta_mat, AclustCPG_df,
                         delta_num, seed_int,
                         betaCols_idx = 9:22,
                         numEx_int = 7, numClusters_int = 500){

  ###  Setup  ###
  # set the seed value from the 'rp'-th element from the vector 'seed_values'
  set.seed(seed_int)

  # randomly pick numClusters_int clusters
  clusts_int <- AclustCPG_df$Clusternumber
  randClusts_int <- sample(max(clusts_int), numClusters_int)

  # full data info of randomly picked clusters
  inflateClust_df <- AclustCPG_df[which(clusts_int %in% randClusts_int), ]
  bval_mat <- as.matrix(inflateClust_df[, betaCols_idx])


  ### Check average b-val/m-val group-wise in each cpg, and add delta value ###
  if (delta_num > 0) {

    # for all cpgs of numClusters_int random picked clusters
    for (i in 1:nrow(bval_mat)) {

      avgbetagr1 <- mean(bval_mat[i, 1:numEx_int])
      avgbetagr2 <- mean(bval_mat[i, (numEx_int + 1):ncol(bval_mat)])

      # Switch statement to ensure that adding delta only increases the signal
      if (avgbetagr1 > avgbetagr2) {
        bval_mat[i, 1:numEx_int] <- bval_mat[i, 1:numEx_int] + delta_num
      } else {

        bval_mat[i, (numEx_int + 1):ncol(bval_mat)] <-
          bval_mat[i, (numEx_int + 1):ncol(bval_mat)] + delta_num

      }

    } # END for()

    # If any entries of the resultant bval_mat, after adding delta_num, are
    #   greater than 1, replace them by 1 (the max beta score)
    bval_mat[bval_mat >= 1] <- 0.999

  } # END if()
  bval_df <- as.data.frame(bval_mat)

  inflateClust_df$actual <- "positive"


  ###  The other clusters  ###
  # Full data info of remaining CPGs, i.e. except randomly picked clusters
  nonInflateClust_df <-
    AclustCPG_df[!(rownames(AclustCPG_df) %in% inflateClust_df$cpg), ]
  nonInflateClust_df$actual <- "negative"

  # Full data info (beta values) of remaining CPGs, i.e. except randomly picked
  #   clusters
  nonInflateBvals_mat <-
    beta_mat[!(rownames(beta_mat) %in% inflateClust_df$cpg), ]


  ###  Combine Treated and Untreated Clusters  ###
  # Combine the CPGs (belonging to selected clusters) and remianing cpgs
  #   (belonging to non-selected clusters)
  treatedAclustCPG_df <- rbind(inflateClust_df, nonInflateClust_df)
  # Combine the treated (delta-value added) CPGs (belonging to the clusters
  #   selected at random) and untreated (delta-value not added) CPGs (belonging
  #   to the non-selected clusters)
  betaResults_df <- rbind(bval_df, nonInflateBvals_mat)


  ###  Order Columns and Rows  ###
  # Reorder the 'actual' column to previous of beta values
  treatedAclustCPG_df <- as.data.frame(treatedAclustCPG_df)
  impVars_df <- treatedAclustCPG_df[c("Clusternumber", "cpg", "CHR", "MAPINFO",
                                      "start_position", "end_position",
                                      "coordinate_37", "chromosome", "actual")]
  otherVars_df <- treatedAclustCPG_df[setdiff(names(treatedAclustCPG_df),
                                              names(impVars_df))]
  treatedAclustCPG_df <- cbind(impVars_df, otherVars_df)

  # Reorder the rows (CPGs) according to clusternumber
  treatedAclustCPGordered_df <-
    treatedAclustCPG_df[order(treatedAclustCPG_df$Clusternumber), ]

  if (delta_num == 0) {
    treatedAclustCPGordered_df$actual <- "negative"
  }


  ###  Return  ###
  list(simBetaVals_df = betaResults_df,
       simAclusters_df = treatedAclustCPGordered_df)

}
