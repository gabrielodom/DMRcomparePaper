#' Plot Venn Diagrams of DMR Overlaps
#'
#' @description Given a directory of best-performing results files from one of
#'    the simulation functions (\code{\link{WriteDMRcateResults}},
#'    \code{\link{WriteProbeLassoResults}},
#'    \code{\link{WriteBumphunterResults}}, or results from the \code{Comb-p}
#'    method in \code{Python}), call the \code{\link{BuildOverlaps}} function
#'    to import the raw data files and DMR overlap lists, then plot those Venn
#'    diagrams and save the plots to a PDF.
#'
#' @param bestResultsDir The name of the directory where the method results from
#'    the best-performing parameter settings are stored. For the full design we
#'    have included (\code{delta = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4)} and
#'    \code{seed = c(100, 210, 330, 450, 680)}), this directory should contain
#'    35 \code{.RDS} files per method.
#' @param figFileName The name of the figure
#' @param device Which graphics device should be used to save the figures?
#'    Defaults to \code{\link[grDevices]{pdf}}. Note that if you use a device
#'    other than PDF (\code{\link[grDevices]{jpeg}} for instance), you can only
#'    plot one \code{delta_num} and \code{seed_int} combination per file.
#' @param plotTitle The title of the plot. This argument is passed to the
#'    \code{\link[ChIPpeakAnno]{makeVennDiagram}} function. The default value of
#'    \code{"default"} will make the plot title "Venn Diagram for mu = DELTA,
#'    rep = INDEX OF SEED".
#' @param delta_num A vector of treatment sizes with values corresponding to one
#'    of the simulations with completed results files in the
#'    \code{bestResultsDir} directory.
#' @param seeds_int A vector of random seeds with values corresponding to one of
#'    the simulations with completed results files in the \code{bestResultsDir}
#'    directory.
#' @param totalTest_int Parameter passed to the
#'    \code{\link[ChIPpeakAnno]{makeVennDiagram}} function. This is an interger
#'    value specifying the total number of tests performed to obtain the list
#'    of peaks. It should be much larger than the number of peaks in the largest
#'    peak set.
#' @param CPGs_df An annotation table that indicates locations of CpGs.
#'    This data frame has CPG IDs as the rows with matching chromosome and
#'    location info in the columns. Specifically, the columns are: \code{ILMNID}
#'     - the CPG ID; \code{chr} - the chromosome label; and \code{MAPINFO} -
#'    the chromosome location. An example is given in the \code{cpgLocation_df}
#'    data set. This data set is only necessary if the results directory
#'    contains Comb-p results with the specified \code{delta} and \code{seed}
#'    values. This is passed to the \code{\link{BuildOverlaps}} function.
#' @param min.cpgs The minimum number of CPGs before we consider a result
#'    significant. Defaults to 5. This argument is only required if the results
#'    directory contains Comb-p results with the specified \code{delta} and
#'    \code{seed} values. This is passed to the \code{\link{BuildOverlaps}}
#'    function.
#' @param ... Dots for additional arguments to be passed to the graphics device
#'
#' @return Nothing. A PDF file of plots is created as a side effect.
#'
#' @importFrom ChIPpeakAnno makeVennDiagram
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   PlotOverlaps(
#'     bestResultsDir = "best_cases_results/",
#'     figFileName = "best_cases_results/resultsFigures/testVenn_allDesigns2"
#'   )
#' }
PlotOverlaps <- function(bestResultsDir,
                         figFileName,
                         device = pdf,
                         plotTitle = "default",
                         delta_num = c(0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4),
                         seeds_int = c(100, 210, 330, 450, 680),
                         totalTest_int = 3063,
                         CPGs_df = cpgLocation_df,
                         min.cpgs = 5,
                         ...){

  # Group by seeds
  design_mat <- expand.grid(seeds_int, delta_num)
  cpg_df <- CPGs_df
  minCPGs <- min.cpgs

  overlapsByMethods_ls <- lapply(1:nrow(design_mat), function(i){

    overlaps_ls <- BuildOverlaps(
      bestResultsDir = bestResultsDir,
      delta = design_mat[i, 2], seed = design_mat[i, 1],
      CPGs_df = cpg_df,
      min.cpgs = minCPGs
    )
    overlapAttr <- attributes(overlaps_ls)

    null_idx <- sapply(overlaps_ls, is.null)
    overlaps_ls <- overlaps_ls[!null_idx]

    attr(overlaps_ls, "delta") <- overlapAttr$delta
    attr(overlaps_ls, "repl")  <- overlapAttr$repl
    overlaps_ls

  })

  CreateHue <- function(n){

    hues = seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]

  }

  device(file = figFileName, ...)
  lapply(overlapsByMethods_ls, function(x){

    delta_num <- attr(x, "delta")
    repl_int  <- attr(x, "repl")

    if(plotTitle == "default"){
      plotTitle <- paste0(
        "Venn Diagram for mu = ", delta_num, ", rep = ", repl_int
      )
    }

    makeVennDiagram(
      x, NameOfPeaks = names(x), totalTest = totalTest_int,
      by = "region", fill = CreateHue(length(x)),
      main = plotTitle
    )

  })

  dev.off()

}
