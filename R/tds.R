#' Plot TDS curves
#'
#' Plots TDS curves based on dominance rates, showing chance and significance lines.
#' @name tds.plot
#' @aliases tds.plot
#' @usage tds.plot(X, attributes = NULL, times = NULL, chance = NULL, signif = NULL,
#'  line.col = 1, lty = 1, lwd = 1, las = 0, xlab = "Time (seconds)",
#'  ylab = "Dominance rate", main = "", height = 8, width = 12, box = FALSE, save.as = "")
#' @param X matrix of dominance rates (Attributes in rows, Times in columns).
#' @param attributes a vector of attribute labels, corresponding to the attributes in \code{X}.
#' @param times a vector of times, corresponding to the times in \code{X}.
#' @param chance proportion indicating the chance level, usually \code{1/length(attributes)} or \code{1/(1+length(attributes))}.
#' @param signif significance level associated with the number of observations and \code{chance}.
#' @param main plot title; see \code{\link[graphics]{plot}}
#' @param xlab,ylab Labels for the x and y axes; see \code{\link[graphics]{plot}}
#' @param line.col A vector of colors for lines corresponding to \code{attributes}; see \code{\link[graphics]{par}}
#' @param lty,lwd line type and weight for attributes; see \code{\link[graphics]{par}}
#' @param las numeric in {0,1,2,3}; the style of the axis labels. See: \code{\link[graphics]{par}}
#' @param height Window height.
#' @param width Window width.
#' @param box draw box around plot area; see: \code{\link[graphics]{box}}
#' @param save.as Filename if the file will be saved.
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # example using 'bars' data set
#' bars.m <- aggregate(bars[, -c(1:4)], list(samples = bars$sample, attribute = bars$attribute), mean)
#' bars.m <- bars.m[order(bars.m$sample, bars.m$attribute), ]
#' attributes <- unique(bars$attribute)
#' times <- get.times(colnames(bars.m)[-c(1:2)])
#' chance <- get.chance(attributes)
#' signif <- get.significance(chance, nrow(unique(bars[, 1:2])))
#' tds.plot(get.smooth(bars.m[bars.m$sample == 1, -c(1:2)]), attributes = attributes,
#'          times = times, chance = chance, signif = signif,
#'          lwd = 2, main = "Bar 1")
#'
#' # it is possible to hide the portion of the plot below the significance line:
#' rect(-2, -0.2, times[length(times)]+2, signif, col = "white", border = "transparent")
#' # re-add axes & significance line
#' axis(1, labels = seq(0, 45, by = 5), at = seq(0, 45, by = 5))
#' axis(2)
#' abline(h=signif, lty=3, col = "grey")
tds.plot <- function(X, attributes = NULL, times = NULL, chance = NULL,
                     signif = NULL, line.col = 1, lty = 1, lwd = 1, las = 0,
                     xlab = "Time (seconds)", ylab = "Dominance rate",
                     main = "", height = 8, width = 12, box = FALSE, save.as = "") {
  requireNamespace("graphics", quietly = TRUE)
  requireNamespace("grDevices", quietly = TRUE)
  if (any(is.na(times))) {
    times <- 1:ncol(X)
  }
  if (any(is.na(attributes))) {
    attributes <- paste0("attribute", as.character(c(1:nrow(X))))
  }
  first.time <- times[1]
  last.time <- times[length(times)]

  if (length(line.col) < length(attributes) & lty == 1) {
    line.col <- pretty_palette(length(attributes))
  }

  grDevices::dev.new(height = height, width = width)

  if(save.as != "") {
    grDevices::postscript(save.as)
  }
  # basic plot area

  graphics::plot(c(first.time, last.time), c(0, 1), type = "n", xlab = xlab,
       ylab = ylab, las = 1, axes = FALSE, main = main)
  graphics::axis(1)
  graphics::axis(2)
  # chance and significance lines
  x.namesel.offset <- 4
  if (!is.na(chance) & !is.na(signif)) {
    if (chance != signif) {
      x.namesel.offset <- 0
    }
  }
  if (!is.na(signif)) {
    graphics::polygon(x = c(first.time - 9.9, first.time - 9.9, last.time +
                    9.9, last.time + 9.9), y = c(0, signif, signif,
                                                 0), col = "gray95", border = NA)
    graphics::abline(h = signif, lty = 2)
    graphics::text("significance", x = x.namesel.offset, y = signif +
           0.01, adj = c(0, 0), col = "gray66", cex = 0.7,
         font = 3)
  }
  if (!is.na(chance)) {
    graphics::abline(h = chance, lty = 2)
    graphics::text("chance", x = 0, y = chance + 0.01, adj = c(0,
                                                     0), col = "gray66", cex = 0.7, font = 3)
  }
  # attribute lines
  for (att in 1:length(attributes)) {
    graphics::lines(x = times, y = X[att, ], col = line.col[att], lty = lty,
          lwd = lwd)
  }
  graphics::legend(x = "topright", legend = attributes,
         text.col = line.col, ncol = ceiling(length(attributes)/6),
         text.font = 3, bty = "n")

  if (box) graphics::box()

  if(save.as != "") grDevices::dev.off()
}

#' Get TDS dominance rates
#'
#' Get TDS dominance rates.
#' @name get.dominance.rates
#' @aliases get.dominance.rates
#' @usage get.dominance.rates(citations, n)
#' @param citations matrix of dominance counts
#' @param n number of observations (evaluations) per cell
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' x <- rbind(c( 6,  6,  8, 14, 16, 22, 22, 21, 13, 11, 14,  7,  7,  6,  5,  3),
#'            c(14, 24, 31, 36, 37, 39, 44, 48, 51, 55, 48, 40, 30, 20, 10,  5),
#'            c( 7,  8,  9, 15, 17, 21, 21, 20, 21, 22, 18, 17, 17, 20, 20, 20),
#'            c(44, 23, 23, 26,  1,  2,  2,  2,  2,  3,  4,  7, 15, 14, 18, 22),
#'            c(20, 30, 20,  0, 20,  7,  2,  0,  4,  0,  7, 20, 22, 31, 38, 41))
#' colnames(x) <- 0:15
#' get.dominance.rates(x, n = 91)
get.dominance.rates <- function(citations, n) {
  citations <- as.data.frame(citations)
    totals <- totals.vec <- as.data.frame(matrix(rep(n,
                                                     times = ncol(citations)), nrow = 1))
    for (r in 1:nrow(citations)) {
      if (r > 1)
        totals <- as.data.frame(rbind(totals, totals.vec))
    }
    dom.props <- citations/totals
  return(dom.props)
}

#' Get times
#'
#' Convenience function to convert exported time labels, e.g. from character format c('time_0.1s', 'time_0.2s', ...) or related format to numeric format c(0.1, 0.2, ...).
#' @name get.times
#' @aliases get.times
#' @usage get.times(time.char, trim.left = "time_", trim.right = "s")
#' @param time.char vector of characters containing the time
#' @param trim.left string to be trimmed from left
#' @param trim.right string to be trimmed from right
#' @return times vector of times in numeric format
#' @export
#' @encoding UTF-8
#' @details Convenience function for getting times from column headers from common data export formats.
#' @examples
#' get.times(colnames(bars)[-c(1:4)])
#'
#' (sample.colnames <- paste0("X", 0:30))
#' get.times(sample.colnames, trim.left = "X", trim.right = "")
get.times <- function(time.char, trim.left = "time_", trim.right = "s") {
  # convenience function converts time labels in format
  # c('time_0.1s', 'time_0.2s', ...) to numeric c(0.1, 0.2,
  # ...)
  for (i in 1:length(time.char)) {
    time.char[i] <- substr(time.char[i], nchar(trim.left) +
                                 1, nchar(time.char[i]) - nchar(trim.right))
  }
  return(times = as.numeric(time.char))
}

#' TDS chance proportion
#'
#' Obtains the TDS chance proportion based on the number of attributes, as proposed by Pineau et al. (2009; Eq. 1).
#' @name get.chance
#' @aliases get.chance
#' @usage get.chance(attributes = c(), include.stop = FALSE)
#' @param attributes number of attributes used in the TDS ballot.
#' @param include.stop defaut is \code{FALSE}. Default should be kept if time standardization is applied. Optionally, set to \code{TRUE} if analyzing data on the raw timeline.
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # example using 'bars' data set
#' attributes <- unique(bars$attribute)
#' chance <- get.chance(attributes)
#' chance
get.chance <- function(attributes = c(), include.stop = FALSE) {
  if (length(attributes) == 0) {
    return(NA)
  } else {
    len <- length(attributes)
    if (include.stop) len <- len + 1
    return(chance = 1 / len)
  }
}

#' TDS significance proportion
#'
#' Obtains the TDS significance proportion based on the number of observations and chance, as proposed by Pineau et al. (2009; Eq. 1).
#' @name get.significance
#' @aliases get.significance
#' @usage get.significance(chance, n, alpha = 0.05)
#' @param chance chance proportion; see \code{\link[tempR]{get.chance}}.
#' @param n number of observations.
#' @param alpha significance level for binomial test of 2 independent proportions (based on normal approximation; see: Pineau et al., 2009, Eq. 1)
#' @details The TDS significance level proposed by Pineau et al. (2009, Eq. 1) provides a simple and widely used heuristic approach for contextualizing observed dominance rates, but should not be used for statistical inference.
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # example using 'bars' data set
#' attributes <- unique(bars$attribute)
#' chance <- get.chance(attributes)
#' signif <- get.significance(chance, nrow(unique(bars[, 1:2])))
#' signif
get.significance <- function(chance, n, alpha = 0.05) {
  if (is.numeric(chance) & is.numeric(n)) {
    requireNamespace("stats", quietly = TRUE)
    signif <- chance + stats::qnorm(1 - alpha) * ((chance * (1 - chance))/n)^(1/2)
  } else {
    signif <- NA
  }
  return(signif = signif)
}

#' Get vector of difference in dominance rates
#'
#' Get vector of difference in dominance rates
#' @name get.differences
#' @aliases get.differences
#' @usage get.differences(x, y)
#' @param x matrix of dominance indicators for a single product\code{x}attribute (rows = evaluations, columns = times)
#' @param y matrix of dominance indicators for a different product (same attribute)
#' @return out vector of differences in dominance rates
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # example using 'bars' data set
#' bars.m <- aggregate(bars[, -c(1:4)], list(samples = bars$sample, attribute = bars$attribute), mean)
#' bars.m <- bars.m[order(bars.m$sample, bars.m$attribute), ]
#' attributes <- unique(bars$attribute)
#' times <- get.times(colnames(bars.m)[-c(1:2)])
#' bar1 <- bars.m[bars.m$sample == 1 & bars.m$attribute == "Caramelized Flavour", -c(1:2)]
#' bar2 <- bars.m[bars.m$sample == 2 & bars.m$attribute == "Caramelized Flavour", -c(1:2)]
#' b.diff <- get.differences(bar1, bar2)
#' round(b.diff, 3)
#'
#' # toy example
#' x <- data.frame(t10 = c( NA,  0,  0,  0,  1,  1,  0,  0,  1,  0, NA),
#'                 t15 = c(  1,  0,  0,  1,  1,  1,  0,  1,  0,  1,  0),
#'                 t20 = c(  1,  1,  1,  1,  1,  1,  1,  0,  1, NA,  0))
#' y <- data.frame(t10 = c( NA, NA,  0,  0,  1,  1,  0,  0,  0,  0, NA),
#'                 t15 = c(  0,  0,  0,  0,  1,  0,  1,  1,  0,  1,  1),
#'                 t20 = c(  1,  0,  1,  1,  0,  0,  1, NA,  1, NA,  0))
#' get.differences(x, y)
get.differences <- function(x, y) {
  return(out = apply(x, 2, mean, na.rm = TRUE) - apply(y, 2, mean, na.rm = TRUE))
}

#' Get least significant differences for pairwise comparisons
#'
#' Get least significant differences for pairwise comparisons (see Pineau et al., 2009, Eq. 2).
#' @name get.significance.diff
#' @aliases get.significance.diff
#' @usage get.significance.diff(x, y, alpha = 0.05)
#' @param x matrix of dominance data (\code{0}/\code{1}) related to one entity
#' @param y matrix of dominance data (\code{0}/\code{1}) related to another entity
#' @param alpha significance for one-sided test (default \code{0.05})
#' @return out least significant difference (at level \code{alpha}) for dominance differences in matrix
#' @details Calculation of least significant differences for TDS difference curves based on Pineau et al. (2009, Eq. 2). The absolute value of the observed dominance rate for a give attribute*time must exceed the corresponding least significant difference calculated here to be considered significant.
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # toy data example
#' x <- data.frame(t10 = c(rep(NA, 15), rep(0, 50), rep(1, 20)),
#'                 t15 = c(rep(NA,  4), rep(0, 61), rep(1, 20)),
#'                 t20 = c(rep(0, 55), rep(1, 30)))
#' y <- data.frame(t10 = c(rep(NA, 15), rep(0, 50), rep(1, 20)),
#'                 t15 = c(rep(NA,  0), rep(0, 21), rep(1, 64)),
#'                 t20 = c( rep(0, 35), rep(1, 50)))
#' signif.xy <- get.significance.diff(x, y)
#' #compare with observed differences
#' diff.xy <- get.differences(x, y)
#' abs(diff.xy) > signif.xy
#'
#' # real data example - differences between Bar 1 and Bar 2 on the attribute "Grain Flavour"
#' attributes <- unique(bars$attribute)
#' times <- get.times(colnames(bars)[-c(1:4)])
#' bar1 <- bars[bars$sample == 1 & bars$attribute == "Grain Flavour", -c(1:4)]
#' bar2 <- bars[bars$sample == 2 & bars$attribute == "Grain Flavour", -c(1:4)]
#' signif.1vs2 <- get.significance.diff(bar1, bar2)
#' # review observed difference in dominance rates vs. least significant differences
#' diff.1vs2 <- get.differences(bar1, bar2)
#' abs(diff.1vs2) > signif.1vs2
#' # differences between samples start at 1.1s and occur throughout the 45.0 evaluation period
get.significance.diff <- function(x, y, alpha = 0.05) {
  if(dim(x)[[2L]] != dim(y)[[2L]]) return("x and y have different numbers of columns")
  x.count <- apply(x, 2, sum, na.rm = TRUE)
  x.n <- nrow(x) - apply(x, 2, lengthwhichis.na)
  y.count <- apply(y, 2, sum, na.rm = TRUE)
  y.n <- nrow(y) - apply(y, 2, lengthwhichis.na)
  P.moy.t = (x.count + y.count)/(x.n + y.n)
  requireNamespace("stats", quietly = TRUE)
  return(out = stats::qnorm(1 - alpha) * sqrt((1/x.n + 1/y.n) * P.moy.t * (1 - P.moy.t)))
}

#' Plot TDS difference curves
#'
#' Plots TDS difference curves based on differences in dominance counts or dominace rates.
#' @name tds.diff.plot
#' @aliases tds.diff.plot
#' @param X matrix of differences in dominance rates (Attributes in rows, Times in columns).
#' @param attributes a vector of attribute labels, corresponding to the attributes in \code{X}.
#' @param times a vector of times, corresponding to the times in \code{X}.
#' @param xlab,ylab Labels for the x and y axes; see \code{\link[graphics]{plot}}
#' @param line.col A vector of colors for lines corresponding to \code{attributes}; see \code{\link[graphics]{par}}
#' @param lty,lwd line type and weight for attributes; see \code{\link[graphics]{par}}
#' @param main plot title; see \code{\link[graphics]{plot}}
#' @details Currently the differences in dominance rates are always displayed. Suppression of differences in dominances rates within a threshold range is not yet implemented.
#' @export
#' @encoding UTF-8
#' @references Pineau, N., Schlich, P., Cordelle, S., Mathonnière, C., Issanchou, S., Imbert, A., Rogeaux, M., Etiévant, P., & Köster, E. (2009). Temporal dominance of sensations: Construction of the TDS curves and comparison with time–intensity.  \emph{Food Quality and Preference}, 20, 450–455. \url{http://dx.doi.org/10.1016/j.foodqual.2009.04.005}
#' @examples
#' # example using 'bars' data set
#' bars.m <- aggregate(bars[, -c(1:4)], list(samples = bars$sample, attribute = bars$attribute), mean)
#' bars.m <- bars.m[order(bars.m$sample, bars.m$attribute), ]
#' attributes <- unique(bars$attribute)
#' times <- get.times(colnames(bars.m)[-c(1:2)])
#' bar1 <- bars.m[bars.m$sample == 1, -c(1:2)]
#' bar2 <- bars.m[bars.m$sample == 2, -c(1:2)]
#' diff.1vs2 <- get.smooth(bar1 - bar2, low.bound = -1, up.bound = 1)
#' tds.diff.plot(diff.1vs2, times = times, attributes = attributes,
#'                 lwd = 2, main = "TDS Differences (Bar 1 - Bar 2)")
#'
#' # suppose we only want to show the curves where the difference in dominance rate
#' # is significantly different
#' # get samples sizes and dominance counts for each product
#' bars.s <- aggregate(bars[, -c(1:4)], list(samples = bars$sample, attribute = bars$attribute), sum)
#' bars.s <- bars.s[order(bars.s$sample, bars.s$attribute), ]
#' bar1.s <- bars.s[bars.s$sample == 1, -c(1:2)]
#' bar2.s <- bars.s[bars.s$sample == 2, -c(1:2)]
#' bar1.n <- nrow(unique(bars[bars$sample == 1, 1:2]))
#' bar2.n <- nrow(unique(bars[bars$sample == 2, 1:2]))
#'
#' # prop.test2 is a wrapper for prop.test (with its default parameters)
#' # thus it will run chi-squared test with Yates continuity correction
#' prop.test2 <-  function(x1, x2, n1, n2, alpha = 0.05){
#'   return((suppressWarnings(prop.test(c(x1,x2), c(n1, n2),
#'           alternative = "two.sided"))$p.value < alpha) %in% TRUE)
#' }
#' # find significant pairwise comparison, treating data as if independent
#' diff_1v2.signif <- mapply(prop.test2, unlist(bar1.s), unlist(bar2.s), bar1.n, bar2.n)
#' # update smoothed difference matrix with NA where non-significant pairs
#' diff_1v2.signif[!diff_1v2.signif] <- NA
#' diff.1vs2 <- diff.1vs2 + diff_1v2.signif - 1
#' # line segments that are non-significant are missing/NA so not plotted
#' tds.diff.plot(diff.1vs2, times = times, attributes = attributes,
#'               lwd = 2, main = "TDS Differences (Bar 1 - Bar 2)")
tds.diff.plot <- function(X, times = NULL, attributes = NULL,
                          xlab = "Time (seconds)", ylab = "Dominance rate",
                          line.col = 1, lty = 1, lwd = 1, main = "") {
  requireNamespace("graphics", quietly = TRUE)
  if (any(is.na(times))) {
    times <- 1:ncol(X)
  }
  if (any(is.na(attributes))) {
    attributes <- paste0("attribute", as.character(c(1:nrow(X))))
  }
  first.time <- times[1]
  last.time <- times[length(times)]
  if (length(line.col) < length(attributes) & lty == 1) {
    line.col <- pretty_palette(length(attributes))
  }
  # basic plot area
  graphics::plot(c(first.time, last.time), c(-1, 1), type = "n", xlab = xlab,
       ylab = ylab, las = 1, main = main)

  # attribute lines
  for (att in 1:length(attributes)) {
    graphics::lines(x = times, y = X[att, ], col = line.col[att], lty = lty,
          lwd = lwd)  # in tds function
  }
  graphics::box()
  graphics::legend(x = "topright", legend = attributes, lwd = lwd, lty = lty,
         col = line.col, ncol = ceiling(length(attributes)/6), cex = 0.7,
         text.font = 1)
}

#' Count observations with missing data
#'
#' Count observations with missing data.
#' @name lengthwhichis.na
#' @aliases lengthwhichis.na
#' @param x  vector data which may contain missings
#' @return \code{count} of observations where data are missing
#' @export
#' @encoding UTF-8
#' @examples
#' x <- c(rep(NA,18), rep(1,18), rep(0,10), rep(NA, 10))
#' lengthwhichis.na(x)
lengthwhichis.na <- function(x){
  return(count = length(which(is.na(x))))
}

#' Time standardize results
#'
#' Set results for a temporal evaluation to a timescale by trimming off time prior to the first onset and following the last offset time, and express the remaining times in terms of percentiles [0, 100].
#' @name std.time
#' @aliases std.time
#' @param X  vector (or data frame) of indicator data.
#' @param trim.left Trim on the left? Default is \code{TRUE}.
#' @param trim.right Trim on the right? Default is \code{TRUE}.
#' @param scale Set to a [0, 1] scale? Default is \code{TRUE}.
#' @param missing indicator for missing data; default is \code{0}.
#' @return out vector (or data frame) of trimmed and/or standardized indicator (\code{0}/\code{1}) data
#' @export
#' @encoding UTF-8
#' @references Lenfant, F., Loret, C., Pineau, N., Hartmann, C., & Martin, N. (2009). Perception of oral food breakdown. The concept of sensory trajectory. \emph{Appetite}, 52, 659-667.
#' @examples
#' # vector - toy data example
#' x <- rep(c(rep(0,18), rep(1,18)), 2)
#' names(x) <- 1:72
#' x           # raw time
#' std.time(x) # standardized time
#'
#' # data frame - toy data example
#' y <- data.frame(rbind(c(c(rep(0,18),
#'                            rep(1,18)),
#'                            rep(0, 4)),
#'                            c(rep(c(rep(0,9),
#'                            rep(1,9)), 2),
#'                            1, rep(0, 3)),
#'                            rep(0, 40)))
#' colnames(y) <- 1:40
#' y           # raw time
#' std.time(y) # standardized time
#'
#' # time standardization using 'bars' data set
#' # only sample 1 will be done (for illustrative purposes)
#' eval1 <- unique(bars[bars$sample == 1, (1:3)])
#' bar1.std <- data.frame(unique(bars[bars$sample == 1, (1:4)]), matrix(0, ncol = 101))
#'
#' for (e in 1:nrow(eval1)){
#'   bar1.std[bar1.std$assessor == eval1$assessor[e] &
#'              bar1.std$session == eval1$session[e] &
#'              bar1.std$sample == eval1$sample[e],
#'              -c(1:4)] <- std.time(bars[bars$assessor == eval1$assessor[e] &
#'                                          bars$session == eval1$session[e] &
#'                                          bars$sample == eval1$sample[e],
#'                                            -c(1:4)])
#' }
#' colnames(bar1.std)[5:ncol(bar1.std)] <- 0:100
#' head(bar1.std)
std.time <- function(X, trim.left = TRUE, trim.right = TRUE, scale = TRUE, missing = 0) {
  if(is.vector(X)){
    it <- 1
    if (all(is.na(X))) return (rep(NA, 101))
    out <- matrix(X, nrow = it)
  } else {
    it <- nrow(X)
    out <- X
    if (all(is.na(X))) return (data.frame(matrix(NA, nrow = it, ncol = 101)))
    if (sum(X) == 0) return (data.frame(matrix(0, nrow = it, ncol = 101)))
  }
  col.y <- ncol(X)

  if (is.na(missing)) {
    out.tot <- out.tot.na <- apply(out, 2, lengthwhichis.na)
    out.tot[out.tot.na == it] <- NA
    out.tot[out.tot.na != it] <- colSums(out[ , out.tot.na != it], na.rm = TRUE)
    ind.keep.start <- min(which(out.tot.na != it))
    ind.keep.stop <- max(which(out.tot.na != it))
  } else {
    out.tot <- colSums(out)
    ind.keep.start <- min(which(out.tot != missing))
    ind.keep.stop <- max(which(out.tot != missing))
  }
  out <- as.matrix(out[ , ind.keep.start:ind.keep.stop], ncol = length(ind.keep.start:ind.keep.stop))
  if (it == 1) out <- matrix(out, nrow = it)
  col.out <- ncol(out)
  if (scale) {
    scale.out <- data.frame(matrix(0, nrow = it, ncol = 101))
    colnames(scale.out) <- 0:100
    rownames(scale.out) <- rownames(X)
    for (i in 1:101) {
      if (i == 1){
        scale.out[, 1] <- out[ , 1]
      } else {
        scale.out[, i] <- out[ , max(1, round(i * col.out/10100, 2)*100)]
      }
    }
    out <- scale.out
  }
  if (it == 1) out <- unlist(c(out))
  return(out)
}

