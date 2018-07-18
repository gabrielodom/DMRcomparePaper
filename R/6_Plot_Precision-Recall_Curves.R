#' Plot Precision-Recall Curves
#'
#' @description Given a list of PR-curve objects as returned by the
#'    \code{\link{BuildPRcurve}} function, plot the precision-recall curve for
#'    each method in a shared figure.
#'
#' @param prCurves_ls A list of PR-curve objects
#' @param new Should the PR curves from this list form their own graph
#'    (\code{TRUE}) or be added onto a previous PR-curve figure (\code{FALSE}).
#'    Defaults to \code{TRUE}.
#' @param lineWidth The line width of each PR curve in the plot. Defaults to 1.
#' @param colours Optionally add your own colours for each line. Otherwise, the
#'    colours are created with the \code{\link[grDevices]{hcl}} function.
#'
#' @return Nothing. A plot is created as a side effect.
#'
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom grDevices hcl
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    prCurves_0.4_100_ls <-
#'      BuildPRcurve(
#'        bestResultsDir = "best_cases_results/",
#'        delta = 0.4,
#'        seed = 100
#'      )
#'
#'    PlotPRCurve(prCurves_0.4_100_ls)
#' }
PlotPRCurve <- function(prCurves_ls,
                        new = TRUE,
                        lineWidth = 1,
                        colours = NULL){

  # Extract Meta
  delta <- attr(prCurves_ls, "delta")
  repl  <- attr(prCurves_ls, "repl")

  # Colours
  if(is.null(colours)){

    CreateHue <- function(n){

      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]

    }

    colours <- CreateHue(length(prCurves_ls))

  }

  # Foundation Plot
  if(new){

    # Main
    plot(
      prCurves_ls[[1]],
      color = colours[1],
      auc.main = FALSE, legend = TRUE,
      main = paste0("Precision-recall curve: mu = ", delta, ", rep = ", repl),
      cex.main = 1, lwd = lineWidth
    )

    # Legend
    legend(x = 0, y = 0.3,
           legend = names(prCurves_ls),
           col = colours, lty = 1, cex = 0.75, box.lty = 0)

  } else {

    plot(
      prCurves_ls[[1]],
      color = colours[1],
      add = TRUE, legend = TRUE,
      lwd = lineWidth
    )

  }


  # Subsequent
  for(i in 2:length(prCurves_ls)){

    if(!is.null(prCurves_ls[[i]])){
      plot(
        prCurves_ls[[i]],
        color = colours[i],
        add = TRUE, legend = TRUE,
        lwd = lineWidth
      )
    }

  }

}
