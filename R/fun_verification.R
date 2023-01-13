#' Verification Rank and Histogram
#'
#' @param forecasts A matrix of ensemble forecasts, in which the rows corresponds to locations and times and the columns correspond to the individual ensemble members.
#' @param observations A vector of observations corresponding to the locations and times of the forecasts.
#' @param do.plot Logical value of plot.
#'
#' @references ensembleBMA::verifRankHist
#'
#' @return A vector giving the rank of verifying observations relative to the corresponding ensemble forecasts. The verification rank historgram is plotted.
#' @export
#'
RankHist <- function (forecasts, observations, do.plot=FALSE)
{
  # rank <- apply(cbind(observations, forecasts), 1, function(x) rank(x, ties = "random",
  #                                                                   na.last = FALSE)[1])
  mat <- as.matrix(cbind(observations, forecasts))
  rank <- matrixStats::rowRanks(mat, ties="random")[,1]

  k <- ncol(forecasts)
  if(do.plot){
    hist(rank, breaks = 0:(k + 1), prob = TRUE, xaxt = "n",
         xlab = "", ylab = "", main = "Verification Rank Histogram")
    axis(1, at = seq(0.5, to = k + 0.5, by = 1), labels = 1:(k + 1))
    abline(h = 1/(k + 1), lty = 2)
  }
  invisible(rank)
}
