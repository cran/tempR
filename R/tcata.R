#' Convenience function for curve smoothing
#'
#' Smooth TCATA curves, constraining smooth within \code{low.bound} and \code{up.bound}.
#' @name get.smooth
#' @aliases get.smooth
#' @usage get.smooth(y, w = NULL, spar = 0.5, low.bound = 0, up.bound = 1)
#' @param y the vector of proportions (or counts) to be smoothed. If a data frame is provided then smoothing is conducted on each row.
#' @param w an optional vector of weights; see \code{\link[stats]{smooth.spline}}
#' @param spar smoothing parameter; see \code{\link[stats]{smooth.spline}}
#' @param low.bound lower bound for smoothed proportions
#' @param up.bound upper bound for smoothed proportions
#' @return out smoothed vector (or data frame with smoothed rows)
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @seealso \code{\link[stats]{smooth.spline}}, \code{\link[stats]{predict}}
#' @examples
#' # example using 'syrah' data set
#' low1 <- t(syrah[seq(3, 1026, by = 6), -c(1:4)])
#' colnames(low1) <- 10:180
#' x <- get.smooth(low1)
#' round(x, 3)
get.smooth <- function(y, w = NULL, spar = 0.5, low.bound = 0, up.bound = 1) {
  if(is.vector(y)){
    it <- 1
    out <- matrix(y, nrow=1)
  } else {
    it <- nrow(y)
    out <- y
  }
  requireNamespace("stats", quietly = TRUE)
  for (i in 1:it){
    incl <- !is.na(out[i, ])
    this.smooth <- stats::smooth.spline(x = seq_along(out[i, ])[incl],
                                        y = out[i, ][incl],
                                 w = w, spar = spar)
    out[i, ] <- stats::predict(this.smooth, x = seq_along(out[i, ]))$y
    out[i, ][out[i, ] < low.bound] <- low.bound
    out[i, ][out[i, ] > up.bound] <- up.bound
  }
  if (it == 1) out <- c(out)
  return(out)
}

#' Fills gaps
#'
#' Replace gaps in TDS and TCATA data with replacement responses.
#' @name fill.gaps
#' @aliases fill.gaps
#' @usage fill.gaps(y, subst = 0, repl = 1)
#' @param y  vector (or data frame) of Bernoulli data which may contain gaps
#' @param subst value occurring in a gap (which represents real data outside a gap). Default is \code{0}.
#' @param repl value occurring for a response (used to replace gap values). Default is \code{1}.
#' @return out vector (or data frame) of Bernoulli data with filled gaps
#' @export
#' @encoding UTF-8
#' @examples
#' # vector with gaps: x with NA gaps (e.g. due to attribute cuing)
#' (x <- rep(c(rep(NA, 4), rep(1, 4)), 2))
#' fill.gaps(x, subst = NA)
#'
#' # array with gaps: y with an gap of 0s (e.g. due to attribute fading)
#' (y <- structure(c(0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0,
#'                   1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0),
#'                 .Dim = c(3L, 10L),
#'                 .Dimnames = list(1:3, 1:10)))
#' fill.gaps(y)
fill.gaps <- function(y, subst = 0, repl = 1) {
  if (all(is.na(c(subst, repl)))) return(y)
  if(is.vector(y)){
    it <- 1
    out <- matrix(y, nrow=1)
  } else {
    it <- nrow(y)
    out <- y
  }
  i.repl.first <- i.repl.last <- i.subst <- NA
  for (i in 1:it){
    if (is.na(subst)){
      i.subst <- which(is.na(out[i, ]))
      ind <- which(out[i, ] == repl)
      if (length(i.subst) > 0 & length(ind) > 0){
        i.repl.first <- min(ind)
        i.repl.last <- max(ind)
        i.subst <- i.subst[i.subst > i.repl.first & i.subst < i.repl.last]
        out[i,][i.subst] <- repl
      }
    } else if (is.na(repl)){
      i.subst <- which(out[i, ] == subst)
      ind <- which(is.na(out[i, ]))
      if (length(i.subst) > 0 & length(ind) > 0){
        i.repl.first <- min(ind)
        i.repl.last <- max(ind)
        i.subst <- i.subst[i.subst > i.repl.first & i.subst < i.repl.last]
        out[i,][i.subst] <- repl
      }
    }  else if (!(is.na(subst) | is.na(subst))){
      if(repl != subst){
        if(any(out[i, ] == repl) & any(out[i, ] == subst)){
          i.repl.first <- min(which(out[i, ] == repl))
          i.repl.last <- max(which(out[i, ] == repl))
          i.subst <- which(out[i, ] == subst)
          i.subst <- i.subst[i.subst > i.repl.first & i.subst < i.repl.last]
          out[i,][i.subst] <- repl
        }
      }
    } else {
      return(warning(print("Could not fill gaps")))
    }
  }
  if (it == 1) out <- c(out)
  return(out)
}

#' Adjust color brightness
#'
#' Select suitable colors for highlighting plots.
#' @name adjust.brightness
#' @aliases adjust.brightness
#' @usage adjust.brightness(rgb.in, percent = 10)
#' @param rgb.in  \code{rgb} of input color
#' @param percent the degree to which input color will be modified/brightened
#' @return hex hex code for new color
#' @export
#' @encoding UTF-8
#' @examples
#' (rgb.in <- c(col2rgb("red")))
#' adjust.brightness(rgb.in, percent = 10)
adjust.brightness <- function(rgb.in, percent = 10) {
  requireNamespace("grDevices", quietly = TRUE)
  col.max <- 253
  col.min <- 0
  new.red <- floor(rgb.in[1] + (100 + percent) / 100)
  new.green <- round(rgb.in[2] + (100 + percent) / 100)
  new.blue <- round(rgb.in[3] + (100 + percent) / 100)
  new.max <- max(new.red, new.green, new.blue)
  new.min <- min(new.red, new.green, new.blue)
  if(new.max >= col.max) {
    # find the ones that are too bright
    col.over <- 1 * (c(new.red, new.green, new.blue) - rep(col.max, times = 3) > 0)
    col.redist <- sum(col.over * (c(new.red, new.green, new.blue) - c(col.max, col.max, col.max)))
    new.red <- ifelse(col.over[1], col.max, min(col.max, round(new.red + col.redist / sum(col.over))))
    # if(col.over[1]){
    #   new.red <- col.max
    # } else {
    #   new.red <- min(col.max, round(new.red + col.redist / sum(col.over)))
    # }
    new.green <- ifelse(col.over[2], col.max, min(col.max, round(new.green + col.redist / sum(col.over))))
    # if(col.over[2]){
    #   new.green <- col.max
    # } else {
    #   new.green <- min(col.max, round(new.green + col.redist / sum(col.over)))
    # }
    new.blue <- ifelse(col.over[3], col.max, min(col.max, round(new.blue + col.redist / sum(col.over))))
    # if(col.over[3]){
    #   new.blue <- col.max
    # } else {
    #   new.blue <- min(col.max, round(new.blue + col.redist / sum(col.over)))
    # }
  }
  if(new.min <= col.min) {
    # one or more colours is too dark
    col.under <- 1 * (c(new.red, new.green, new.blue) - c(col.min, col.min, col.min) < 0)
    col.redist <- sum(col.under * (rep(col.min, times = 3) - c(new.red, new.green, new.blue)))
    new.red <- ifelse(col.under[1], col.min, max(col.min, round(new.red + col.redist / sum(col.under))))
    # if(col.under[1]){
    #   new.red <- col.min
    # } else {
    #   new.red <- max(col.min, round(new.red + col.redist / sum(col.under)))
    # }
    new.green <- ifelse(col.under[2],
                        col.min,
                        max(col.min, round(new.green + col.redist /
                                             sum(col.under))))
    # if(col.under[2]){
    #   new.green <- col.min
    # } else {
    #   new.green <- max(col.min, round(new.green + col.redist / sum(col.under)))
    # }
    new.blue <- ifelse(col.under[3], col.min,
                       max(col.min,
                           round(new.blue + col.redist / sum(col.under))))
    # if(col.under[3]){
    #   new.blue <- col.min
    # } else {
    #   new.blue <- max(col.min, round(new.blue + col.redist / sum(col.under)))
    # }
  }
  return(grDevices::rgb(new.red, new.green, new.blue, max = 255))
}

#' Get decluttering matrix indicating where to show/hide reference lines
#'
#' Declutter TCATA curves by hiding reference lines from plots showing TCATA curves.
#' @name get.decluttered
#' @aliases get.decluttered
#' @usage get.decluttered(x = x, n.x = n.x, y = y, n.y = n.y, alpha = 0.05)
#' @param x selections for sample of interest (can be a vector if several samples of interest)
#' @param n.x evaluations of \code{x} (can be a vector if several samples of interest)
#' @param y selections for comparison (can be a vector if several comparisons will be made)
#' @param n.y evaluations of \code{y} (can be a vector if several comparisons of interest)
#' @param alpha significance level
#' @return declutter vector in which \code{1} indicates "show" and \code{NA} indicates "hide"
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @seealso \code{\link[stats]{fisher.test}}, \code{\link[tempR]{citation.counts}}
#' @examples
#'
#' # functionality of get.decluttered() is conveniently provided in citation.counts()
#'
#' # Data set: ojtcata
#' # Get declutter matrix for comparison of Product 2 vs. average of all products
#' data(ojtcata)
#' oj2.v.all <- citation.counts(ojtcata, product.name = "2", product.col = 2,
#'        attribute.col = 4, results.col = 5:25, comparison = "average")
#' oj2.v.all$declutter
#'
#' # same as
#'
#' p2.declutter <- get.decluttered(x = c(oj2.v.all$P1), n.x = oj2.v.all$Pn,
#'                                 y = c(oj2.v.all$C1), n.y = oj2.v.all$Cn)
#' (p2.declutter <- matrix(p2.declutter, nrow = nrow(oj2.v.all$P1)))
get.decluttered <- function(x = x, n.x = n.x, y = y, n.y = n.y, alpha = 0.05){
  if(any(is.na(x), is.na(y), is.na(n.x), is.na(n.y))) return(print("All parameters required"))
  if(length(x) != length(y)) return(print("Length of x and y must be equal"))
  if(length(n.x) > 1 & length(x) != length(n.x)) return(print("Length of x and n.x are mismatched"))
  if(length(n.y) > 1 & length(y) != length(n.y)) return(print("Length of y and n.y are mismatched"))
  tmp <- data.frame(x1 = x, x0 = n.x - x, y1 = y, y0 = n.y - y)
  fisher.test2 <- function(x){
    requireNamespace("stats", quietly = TRUE)
    return(stats::fisher.test(matrix(x, nrow = 2))$p)
  }
  declutter =  1 * (apply(tmp, 1, fisher.test2) < alpha)
  declutter[declutter == 0] <- NA
  return(declutter = declutter)
}

#' Counts TCATA Citations and Observations for a Product and a Comparison Set
#'
#' Calculates how many times a specified product was checked and how many times a comparison set was checked.
#' The number of evaluations for the product and comparison set are also calculated,
#' along with a reference and decluttering matrix for plotting in \code{\link[tempR]{tcata.line.plot}}.
#' @name citation.counts
#' @aliases citation.counts
#' @usage citation.counts(x, product.name = "", product.col = 1,
#' attribute.col = 2, results.col = NULL, comparison = "average")
#' @param x matrix of TCATA 0/1 data with (Assessors x Products x Reps x Attributes) in rows with row headers and (Times) in columns
#' @param product.name name of the product for which to calculate how many times a product was checked and not checked
#' @param product.col index of column in \code{x} that contains the product identities
#' @param attribute.col index of column in \code{x} that contains the attribute identities
#' @param results.col indices of columns in \code{x} that contain the raw (0/1) data
#' @param comparison specifies whether the comparison will be with the average of \emph{all} products (\code{"average"} (default)) or with the average of the \emph{other} products (\code{"other"}, i.e. excludes the product specified by \code{product.name})
#' @return list object with three elements:
#' \itemize{
#' \item{\code{P1} matrix of counts for product specified by \code{product.name}
#' (attributes are in rows; times are in columns).}
#' \item{\code{Pn} number of observations for \code{product.name}}
#' \item{\code{C1} matrix of counts for comparison set specified by \code{comparison}
#' (dimensions equal to \code{P1}.}
#' \item{\code{Cn} number of observations for the comparison set defined by \code{comparison}}
#' \item{\code{ref} a matrix of citation proportions for the comparison set specified
#' by \code{comparison} (dimensions equal to \code{P1}; can be used to draw a reference line;
#' see \code{\link[tempR]{tcata.line.plot}}}
#' \item{\code{declutter} a matrix for decluttering in a line plot
#' (dimensions equal to \code{P1}; see \code{\link[tempR]{get.decluttered}}}}
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @references Meyners, M., Castura, J.C. (2018). The analysis of temporal check-all-that-apply (TCATA) data. \emph{Food Quality and Preference}, 67, 67-76. \doi{10.1016/j.foodqual.2017.02.003}
#' @seealso \code{\link[tempR]{tcata.line.plot}}, \code{\link[tempR]{get.decluttered}}
#' @examples
#' # example using 'ojtcata' data set
#' data(ojtcata)
#'
#' # comparison of Orange Juice 3 vs. all other OJs (1, 2, 4, 5, 6)
#' oj3.v.other <- citation.counts(ojtcata, product.name = "3", product.col = 2,
#'        attribute.col = 4, results.col = 5:25, comparison = "other")
#'
#' # show results
#' oj3.v.other
#'
#' times <- get.times(colnames(ojtcata)[-c(1:4)])
#' attributes <- unique(ojtcata$attribute)
#' palettes <- make.palettes(length(attributes))
#'
#' # plot results
#' tcata.line.plot(oj3.v.other$P1, n = oj3.v.other$Pn,
#'    attributes = attributes, times = times,
#'    line.col = palettes$pal, reference = oj3.v.other$ref, ref.lty = 3,
#'    declutter = oj3.v.other$declutter, highlight = TRUE, highlight.lwd = 4,
#'    highlight.col = palettes$pal.light,
#'    height = 7, width = 11, legend.cex = 0.7, main = "Product 3 vs. Other Products")
citation.counts <- function(x, product.name = "",
                          product.col = 1, attribute.col = 2,
                          results.col = NULL, comparison = "average"){
  if(product.name == ""){
    product.name <- unique(x[, product.col])[1]
  }
  if(is.null(results.col)){
    results.col <- (1:ncol(x))[-c(product.col, attribute.col)]
  }
  AllProducts1 <- stats::aggregate(x[, results.col],
                            list(samp = x[, product.col], attribute = x[, attribute.col]),
                            sum)
  P1 <- as.matrix(AllProducts1[AllProducts1$samp == product.name, -c(1:2)])
  Pn <- nrow(x[x[,product.col] == product.name &
                 x[,attribute.col] == unique(x[,attribute.col])[1], ])
  if(comparison == "average"){
    C1 <- as.matrix(stats::aggregate(AllProducts1[, -c(1:2)],
                    list(attribute = AllProducts1[, 2]),
                    sum)[,-1])
    Cn <- nrow(x)/length(unique(x[,attribute.col]))
  }
  if(comparison == "other"){
    C1 <- as.matrix(stats::aggregate(AllProducts1[AllProducts1$samp !=
                                                    product.name, -c(1:2)],
                    list(attribute = AllProducts1[AllProducts1$samp !=
                                                    product.name, 2]),
                    sum)[,-1])
    Cn <- nrow(x[x[,product.col] != product.name &
                   x[,attribute.col] == unique(x[,attribute.col])[1], ])
  }
  ref <- C1 / Cn
  declutter <-  matrix(get.decluttered(x = c(P1), n.x = Pn,
                                       y = c(C1), n.y = Cn),
                       nrow = nrow(P1), dimnames = dimnames(P1))

  return(list(P1 = P1, Pn = Pn, C1 = C1, Cn = Cn,
              ref = ref, declutter = declutter))
}

#' Temporal Check-All-That-Apply (TCATA) curve
#'
#' Plots TCATA curves based on count or proportion data. Can also be used for plotting Temporal Dominance of Sensations (TDS) curves based on dominance counts or proportions.
#' @name tcata.line.plot
#' @aliases tcata.line.plot
#' @usage tcata.line.plot(X, n = 1, attributes = c(), times = c(),
#' lwd = 1, lty = 1, line.col = c(),
#' emphasis = NA, emphasis.col = c(), emphasis.lty = 1, emphasis.lwd = 3,
#' declutter = NA,
#' reference = NA, ref.col = c(), ref.lty = 2, ref.lwd = 1,
#' highlight = FALSE, highlight.col = c(), highlight.lty = 1, highlight.lwd = 5,
#' xlab = "Time", ylab = "Citation proportion", axes.font = 1,
#' axes.cex = 1, xlim = c(), las = 0,
#' x.increment = 5, box = FALSE,
#' legend.cex = 1, legend.font = 1, legend.pos = "topleft", legend.ncol = 2,
#' height = 8, width = 12, main = "",
#' save.format = "", save.as = "" )
#' @param X matrix of proportions (or, if there is no missing data, on counts), typically with Attributes in rows and times in columns.
#' @param n The number of observations if \code{X} is a count matrix. Keep \code{n = 1} if \code{X} is a matrix of proportions.
#' @param attributes a vector of attribute labels, corresponding to the attributes in \code{X}.
#' @param times a vector of time, corresponding to the times in \code{X}.
#' @param lwd line width for attribute curves that matches either \code{attributes} or \code{X}.
#' @param lty line types for attribute curves that matches either \code{attributes} or \code{X}.
#' @param line.col attribute curves colours that matches \code{attributes}.
#' @param emphasis matrix matching \code{X} in its dimensions, with a numeric value corresponding to points requiring emphasis, and \code{NA} for points without emphasis.
#' @param emphasis.col vector colours for attributes corresponding to rows of \code{X}; taken from \code{line.col} if not specified.
#' @param emphasis.lty either a line type (\code{lty}) for all emphasis lines .
#' @param emphasis.lwd line weight associated with the emphasis line.
#' @param declutter a matrix with the same dimensions as \code{X}; give the value \code{1} to show a proportion in \code{X} and \code{reference} (if given), otherwise give \code{0} or \code{NA}.
#' @param reference a matrix with the same dimensions as \code{X}; give the value \code{1} if \code{reference} will be shown (allowing finer control than \code{declutter}), otherwise give \code{0} or \code{NA}
#' @param ref.col \code{reference} line colour
#' @param ref.lty \code{reference} line type
#' @param ref.lwd \code{reference} line width
#' @param highlight TRUE if differences will be highlighted; otherwise FALSE
#' @param highlight.col a vector of colours for attributes corresponding to rows of \code{X}
#' @param highlight.lty line type associated with the highlighting
#' @param highlight.lwd line weight associated with the highlighting line
#' @param xlab label for the x axis
#' @param ylab label for the y axis
#' @param axes.font font for axes labels; see \code{\link[graphics]{par}}
#' @param axes.cex size for axes labels.
#' @param xlim x limits specified using a vector of 2 (ascending) numbers.
#' @param las numeric in \code{0,1,2,3} indicating style of axis labels; see \code{\link[graphics]{par}}
#' @param x.increment interval between times when labelling the x axis
#' @param box draw box around plot area; see: \code{\link[graphics]{box}}
#' @param legend.cex size of markers shown in the legend
#' @param legend.font font for the legend; see \code{\link[graphics]{text}}
#' @param legend.pos location of plot legend; defaults to \code{"topleft"}
#' @param legend.ncol number of columns in legend
#' @param main plot title; see \code{\link[graphics]{plot}}
#' @param height window height
#' @param width window width
#' @param save.format If indicated, this will be the fle type for the save image. Defaults to \code{"eps"} (eps format). Other possible values are \code{""} (not saved) or \code{"png"} (png format)
#' @param save.as Filename if the file will be saved
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @references Meyners, M., Castura, J.C. (2018). The analysis of temporal check-all-that-apply (TCATA) data. \emph{Food Quality and Preference}, 67, 67-76. \doi{10.1016/j.foodqual.2017.02.003}
#' @examples
#' # example using 'syrah' data set
#' low1 <- t(syrah[seq(3, 1026, by = 6), -c(1:4)])
#' colnames(low1) <- 10:180
#' tcata.line.plot(get.smooth(low1), lwd = 2, main = "Low-ethanol wine (Sip 1)")
#'
#' # example using 'ojtcata' data set
#' data(ojtcata)
#' # comparison of Orange Juice 1 vs. Other OJs (2 to 6)
#' oj1.v.other <- citation.counts(ojtcata, product.name = "1", product.col = 2,
#'        attribute.col = 4, results.col = 5:25, comparison = "other")
#' times <- get.times(colnames(ojtcata)[-c(1:4)])
#' attributes <- unique(ojtcata$attribute)
#' palettes <- make.palettes(length(attributes))
#'
#' # plot results
#' tcata.line.plot(oj1.v.other$P1, n = oj1.v.other$Pn,
#'    attributes = attributes, times = times,
#'    line.col = palettes$pal, reference = oj1.v.other$ref, ref.lty = 3,
#'    declutter = oj1.v.other$declutter, highlight = TRUE, highlight.lwd = 4,
#'    highlight.col = palettes$pal.light,
#'    height = 7, width = 11, legend.cex = 0.7, main = "Product 1 vs. Other Products")
#'
#' # example showing plots similar to those in Meyners & Castura (2018)
#' # comparison of Orange Juice 1 vs. All OJs (1 to 6)
#' oj1.v.all <- citation.counts(ojtcata, product.name = "1", product.col = 2,
#'        attribute.col = 4, results.col = 5:25, comparison = "average")
#' lty.mat <- matrix(1,nrow=6,ncol=21)
#' lty.mat[, 1:3] <- c(rep(NA,8),rep(c(1,NA),4), 1, 1)
#' lty.mat[2, 9:12] <- lty.mat[5, 8] <- 3
#' tcata.line.plot(oj1.v.all$P1, n = oj1.v.all$Pn, attributes = attributes,
#'                 times = times, line.col = palettes$pal, lty = lty.mat, lwd = 2,
#'                 height = 7, width = 11, legend.cex = 0.7, main = "Product 1 vs. All Products")
tcata.line.plot <- function(X, n = 1, attributes = c(), times = c(),
                            lwd = 1, lty = 1, line.col = c(),
                            emphasis = NA, emphasis.col = c(),
                            emphasis.lty = 1, emphasis.lwd = 3,
                            declutter = NA,
                            reference = NA, ref.col = c(),
                            ref.lty = 2, ref.lwd = 1,
                            highlight = FALSE, highlight.col = c(),
                            highlight.lty = 1, highlight.lwd = 5,
                            xlab = "Time", ylab = "Citation proportion",
                            axes.font = 1, axes.cex = 1,
                            xlim = c(), las = 0, x.increment = 5, box = FALSE,
                            legend.cex = 1, legend.font = 1,
                            legend.pos = "topleft", legend.ncol = 2,
                            height = 8, width = 12, main = "",
                            save.format = "eps", save.as = "" ){

  # mat will contain proportions
  requireNamespace("grDevices", quietly = TRUE)
  requireNamespace("graphics", quietly = TRUE)
  if (length(attributes) == 0 ) attributes <- rownames(X)
  if (length(times) == 0 ) times <- as.numeric(colnames(X))
  if (length(line.col) == 0 ) line.col <- pretty_palette(length(attributes))
  if (length(highlight.col) == 0 ){
    highlight.col <- make.palettes(length(attributes))$pal.light
  }

  X <- X/n
  grDevices::dev.new(height = height, width = width)
  start.time <- min(times)
  end.time <- max(times)
  if(length(xlim) == 0){
    xlim.given <- FALSE
    xlim <- c(start.time, end.time)
    xlim.at <- seq(from = start.time, to = end.time, by = x.increment)
  } else {
    xlim.given <- TRUE
    xlim.at <- seq(from = start.time, to = end.time, by = x.increment)
    if(xlim.at[1] != xlim[1]) xlim.at <- c(xlim[1], xlim.at)
    if(xlim.at[length(xlim.at)] != xlim[2]) xlim.at <- c(xlim.at, xlim[2])
  }

  for(gg in 1:length(save.format)){
    if(save.format[gg] == "eps" & save.as[gg] != "") {
      grDevices::postscript(save.as[gg])
    }

    graphics::plot(x = c(start.time, end.time), y = c(0,1),
                   xlim = xlim, ylim = c(0,1),
         xlab = xlab, ylab = ylab, font.lab = axes.font, cex.lab = axes.cex,
         type = "n", axes = FALSE, main = main)
    graphics::axis(1, labels = xlim.at, at = xlim.at)
    graphics::axis(2)

    if (!any(is.na(reference))) {
      for (a in seq_along(attributes)){
        declutter.a <- rep(1, length = length(reference[a, ]))
        if(length(1*is.na(declutter)) != 1){
          declutter.a <- declutter[a, ]
        }
        graphics::lines(x = times, y = reference[a, ]*declutter.a,
                        xlim = c(start.time, end.time),
                        col = line.col[a], lwd = ref.lwd, lty = ref.lty)
        if(highlight){
          graphics::lines(x = times, y = X[a, ]*declutter.a,
                          xlim = c(start.time, end.time),
                          col = highlight.col[a],
                          lwd = highlight.lwd, lty = highlight.lty)
        }
      }
    }

    if(!any(is.na(emphasis))){
      if (!is.na(nrow(emphasis))) {
        if (nrow(emphasis) == length(attributes)){
          if(length(emphasis.col)==0) emphasis.col <- line.col
          for (a in seq_along(attributes)){
            emphasis[a,][emphasis[a,] == 0] <- NA
            graphics::lines(x = times, y = X[a,]*emphasis[a,],
                            xlim = c(start.time, end.time),
                            col = emphasis.col[a], lwd = emphasis.lwd,
                            lty = emphasis.lty)
          }
        }
      }
    }

    if ((length(lwd) == 1) && (length(lty) == 1)){
      for (a in seq_along(attributes)){
        graphics::lines(x = times, y = X[a,], xlim = c(start.time, end.time),
                        col = line.col[a], lwd = lwd, lty = lty)
      }
    } else {
      if (length(lwd) == 1) lwd <- matrix(lwd, nrow = nrow(X), ncol = ncol(X))
      if (length(lty) == 1) lty <- matrix(lty, nrow = nrow(X), ncol = ncol(X))
      for (a in seq_along(attributes)){
        aX <- X[a, ]
        lstyle.unique <- unique(lstyle <- data.frame(lty = c(t(lty[a, ])),
                                                     lwd = c(t(lwd[a, ]))))
        for (this.lstyle in 1:nrow(lstyle.unique)){
          declutter.a <- rep(1, length = length(aX))
          if (!any(is.na(reference))) {
            if(length(1*is.na(declutter)) != 1){
              declutter.a <- declutter[a, ]
            }
          }
          ko.X <- data.frame(X = unlist(aX), lty = lstyle$lty,
                             lwd = lstyle$lwd, declutter = declutter.a, ko = 1)
          # draw the line between the 1st item in the list and the last
          ko.X$ko[is.na(declutter.a) & is.na(c(declutter.a[-1],1))] <- NA
          ko.X$ko[lstyle$lty != lstyle.unique$lty[this.lstyle]] <- NA
          ko.X$ko[lstyle$lwd != lstyle.unique$lwd[this.lstyle]] <- NA
          ko.X$ko[is.na(lstyle$lty)] <- NA
          ko.X$ko[is.na(lstyle$lwd)] <- NA
          if(!is.na(lstyle.unique$lwd[this.lstyle]) &
             !is.na(lstyle.unique$lty[this.lstyle])){
            graphics::lines(x = times,
                            y = ko.X$X * ko.X$ko,
                            col = line.col[a],
                            lwd = lstyle.unique$lwd[this.lstyle],
                            lty = lstyle.unique$lty[this.lstyle])
          }
        }
      }
    }

    graphics::legend(legend.pos, legend = attributes, bty = "n",
                     ncol = legend.ncol, text.col = line.col,
                     text.font = legend.font, cex = legend.cex)

    if (box) graphics::box()

    if(save.format[gg] == "png" & save.as[gg] != ""){
      grDevices::savePlot(save.as[gg], type = "png")
    }
    if(save.format[gg] == "eps" & save.as[gg] != ""){
      grDevices::dev.off()
    }
  }
}

#' Pairwise comparisons
#'
#' p-value for pairwise comparisons.
#' @name get.mat.diff.sign
#' @aliases get.mat.diff.sign
#' @param x citations for product x
#' @param y citations for product y
#' @param n.x total observations for x
#' @param n.y total observations for y
#' @param test.type So far only Fisher's exact test is implemented (\code{"f"})
#' @seealso \code{\link[stats]{fisher.test}}
#' @aliases get.mat.diff.sign
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @examples
#' # Toy TCATA citations data for two samples: s1, s2
#' s1 <- t(data.frame(sweet =  c(10, 23, 25, 26, 26, 43, 44),
#'                    bitter = c( 4, 18, 19, 27, 36, 43, 54),
#'                    sour = c(40, 53, 85, 70, 46, 33, 24)))
#' s2 <- t(data.frame(sweet = c(11, 33, 45, 46, 56, 43, 44),
#'                    bitter = c( 0, 11, 11, 14, 25, 35, 34),
#'                    sour = c(30, 33, 35, 20, 26, 23, 24)))
#' colnames(s1) <- colnames(s2) <- paste0("time_", seq(5, 35, by = 5), "s")
#' n <- 90
#' signif <- get.mat.diff.sign(s1, s2, n, n)
#' signif
get.mat.diff.sign <- function(x = x, y = y, n.x = n.x, n.y = n.x, test.type = "f"){
  requireNamespace("stats", quietly = TRUE)
  x1 <- unlist(c(x))
  y1 <- unlist(c(y))
  out.mat <- x * 0
  tmp <- data.frame(x1 = x1, x0 = c(n.x) - x1, y1 = y1, y0 = c(n.y) - y1)
  fisher.test2 <- function(X){
    requireNamespace("stats", quietly = TRUE)
    return(stats::fisher.test(matrix(X, nrow = 2))$p)
  }
  out.mat <- matrix(apply(tmp, 1, fisher.test2),
                    nrow = nrow(as.matrix(x)), ncol = ncol(as.matrix(x)))
  return(out.mat)
}

#' TCATA difference plot
#'
#' Plots TCATA difference curves.
#' @name tcata.diff.plot
#' @aliases tcata.diff.plot
#' @usage tcata.diff.plot(x1 = x1, x2 = NA, n1 = 1, n2 = NA,
#' attributes = c(), times = c(), lwd = 1,
#' declutter = NA, get.decluttered = FALSE, emphasis = NA, alpha = 0.05, emphasis.lwd = 3,
#' main = "", height = 8, width = 12,
#' xlab = "Time", ylab = "Difference in citation proportion",
#' axes.font = 1, axes.cex = 1, line.col = c(), x.increment = 5,
#' legend.cex = 1, legend.font = 1, save.as = "")
#' @param x1 matrix of difference proportions, or of counts if \code{n1} specified. If \code{mat2} specified then proportions or counts apply to first sample. Attributes are in rows, times in columns.
#' @param x2 matrix of proportions for second sample, or of counts if \code{n2} specified.
#' @param n1 number of observations for first sample
#' @param n2 number of observations for second sample
#' @param attributes vector of attribute labels for row in \code{x1} (and \code{x2})
#' @param times vector of times for columns in \code{x1} (and \code{x2})
#' @param lwd Line width
#' @param declutter indicator matrix with same dimensions of \code{x1} to suppress output
#' @param get.decluttered if \code{TRUE} then calculates the \code{declutter}  matrix from \code{get.mat.diff.sign}
#' @param emphasis set to \code{1} to emphasize significant differences
#' @param alpha significance level for entrywise test of \code{x1} and \code{x2} (if counts)
#' @param emphasis.lwd line weight for emphasizing significant differences
#' @param main plot title; see \code{\link[graphics]{plot}}
#' @param height plot height
#' @param width plot width
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param axes.font Font for axes labels; see \code{\link[graphics]{par}}.
#' @param axes.cex Size for axes labels.
#' @param line.col line color for attribute lines
#' @param legend.cex symbol size for legend
#' @param legend.font Font for the legend; see \code{\link[graphics]{text}}.
#' @param x.increment increment between time labels on x axis
#' @param save.as Filename to use if file will be saved.
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @examples
#' # difference between High and Low ethanol wines (sip 1)
#' x.diff.raw <- t(syrah[seq(1, 1026, by = 6), -c(1:4)]) -
#'                 t(syrah[seq(3, 1026, by = 6), -c(1:4)])
#' x.diff.smooth <- get.smooth(x.diff.raw, low.bound = -1, up.bound = 1)
#' colnames(x.diff.smooth) <- colnames(x.diff.raw) <- times <- 10:180
#' tcata.diff.plot(x1 = x.diff.smooth, attributes = rownames(x.diff.smooth), times = times, lwd = 2,
#'                 main = "Sip 1 differences: High-ethanol wine - Low-ethanol wine")
#'
#' # an example based on the syrah data set (truncated for efficiency)
#' n <- 52
#' H1 <- t(syrah[seq(1, 126, by = 6), -c(1:4)] * n)
#' L1 <- t(syrah[seq(3, 126, by = 6), -c(1:4)] * n)
#' colnames(H1) <- colnames(L1) <- times <- 10:30
#' tcata.diff.plot(x1 = H1, x2 = L1, n1 = n, n2 = n,
#'                 attributes = rownames(H1), get.decluttered = TRUE, lwd = 2)
tcata.diff.plot <- function(x1 = x1, x2 = NA, n1 = 1, n2 = NA, attributes = c(), times = c(), lwd = 1,
                            declutter = NA, get.decluttered = FALSE, emphasis = NA, alpha = 0.05, emphasis.lwd = 3,
                            main = "", height = 8, width = 12,
                            xlab = "Time", ylab = "Difference in citation proportion",
                            axes.font = 1, axes.cex = 1, line.col = c(),
                            x.increment = 5, legend.cex = 1, legend.font = 1,
                            save.as = ""){
  requireNamespace("grDevices", quietly = TRUE)
  requireNamespace("graphics", quietly = TRUE)
  # expect x1 and x2 to be contigency tables with the same dimension, rownames, colnames
  if (sum(1*(dim(x1) != dim(x2))) + sum(1*(colnames(x1) != colnames(x2))) != 0){
    print("Mismatch in matrix dimension or matrix attributes (colnames).")
    return(NULL)
  }
  if (length(times) == 0) times <- as.numeric(colnames(x1))
  if (length(attributes) == 0) attributes <- rownames(x1)
  if (length(line.col) == 0) {
    requireNamespace("grDevices", quietly = TRUE)
    line.col <- pretty_palette(length(attributes))
  }
  mat.diff.sign <- NA
  if (!all(is.na(c(x2, n2)))){
    # two matrices are provided - get the differences
    mat.diff <- (x1/n1) - (x2/n2)
    if(get.decluttered){
      # override any declutter matrix that was supplied
      mat.diff.sign <- get.mat.diff.sign(x1, x2, n1, n2)
      mat.diff.sign[mat.diff.sign > alpha] <- NA # hide values that are nsd
      mat.diff.sign[!is.na(mat.diff.sign)] <- 1
      colnames(mat.diff.sign) <- colnames(mat.diff)
    }
  } else {
    # only one matrix of differences
    mat.diff <- x1 / n1
  }

  grDevices::dev.new(height = height, width = width)
  start.time <- min(times)
  end.time <- max(times)

  if(save.as != "") {
    grDevices::postscript(save.as)
  }

  graphics::plot(x = c(start.time, end.time), y = c(-1, 1), axes = FALSE, xlab = xlab, ylab = ylab,
                 font.lab = axes.font, cex.lab = axes.cex, type = "n", main = main)
  graphics::axis(1, at = seq(from = min(as.numeric(colnames(x1))), to = max(as.numeric(colnames(x1))), by = x.increment))
  graphics::axis(2)
  graphics::abline(h = 0, col = "grey33", lty = 3)
  if(all(dim(x1) == dim(x2)) & !(all(is.na(x2)))){
    if((get.decluttered || !all(is.na(mat.diff.sign)))){
      if(all(dim(mat.diff) == dim(mat.diff.sign))){
        # declutter the plot
        for (a in seq_along(attributes)){
          for (t in 1:ncol(mat.diff.sign)){
            continue <- FALSE
            start.x <- stop.x <- NA
            if(!is.na(mat.diff.sign[a, t])){
              if(t == 1) continue <- TRUE
              if(continue==FALSE){
                if(is.na(mat.diff.sign[a, t - 1])) continue <- TRUE
              }
              if(continue){
                start.x <- t
                findnext0 <- which(is.na(mat.diff.sign[a,min(start.x+1,length(colnames(mat.diff.sign))):length(colnames(mat.diff.sign))]))
                # use the first time that is non-significant, or the last time
                stop.x <- ifelse(sum(findnext0*1)==0, length(colnames(mat.diff.sign)), findnext0[1] - 1 + start.x)
                #if (sum(findnext0*1)==0) {
                #  stop.x <- length(colnames(mat.diff.sign.out)) # all significant;
                #} else {
                #  stop.x <- findnext0[1] - 1 + start.x # take the first time that is non-significant
                #}
                y.start.x <- start.x
                y.stop.x <- stop.x
                graphics::lines(x = as.numeric(colnames(mat.diff)[start.x:stop.x]), y = mat.diff[a, c(y.start.x:y.stop.x)], col = line.col[a], lwd = lwd)
              }
            }
          }
        }
      }
    }
  } else {
    for (a in seq_along(attributes)){
      graphics::lines(x = times, y = mat.diff[a,], lwd = lwd, col = line.col[a])
    }
  }
  if(!any(is.na(emphasis))){
    if (!is.na(nrow(emphasis))) {
      if (nrow(emphasis) == length(attributes)){
        # emphasize significant differences
        mat.diff.sign <- get.mat.diff.sign(x1, x2, n1, n2)
        mat.diff.sign[mat.diff.sign > alpha] <- NA # hide values that are nsd
        colnames(mat.diff.sign) <- colnames(mat.diff)
        mat.diff.sign.out <- mat.diff.sign

        for (a in seq_along(attributes)){
          for (t in seq_along(colnames(mat.diff.sign.out))){
            continue <- FALSE
            start.x <- stop.x <- NA
            if(!is.na(mat.diff.sign.out[a, t])){
              if(t == 1) continue <- TRUE
              if(continue==FALSE){
                if(is.na(mat.diff.sign.out[a, t - 1])) continue <- TRUE
              }
              if(continue){
                start.x <- t
                findnext0 <- which(is.na(mat.diff.sign.out[a,min(start.x+1,length(colnames(mat.diff.sign.out))):length(colnames(mat.diff.sign.out))]))
                # use the first time that is non-significant, or the last time
                stop.x <- ifelse(sum(findnext0*1)==0, length(colnames(mat.diff.sign.out)), findnext0[1] - 1 + start.x)
                #if (sum(findnext0*1)==0) {
                #  stop.x <- length(colnames(mat.diff.sign.out)) # all significant;
                #} else {
                #  stop.x <- findnext0[1] - 1 + start.x # take the first time that is non-significant
                #}
                y.start.x <- start.x
                y.stop.x <- stop.x
                graphics::lines(x = as.numeric(colnames(mat.diff)[start.x:stop.x]), y = mat.diff[a, c(y.start.x:y.stop.x)], col = line.col[a], lwd = emphasis.lwd)
              }
            }
          }
        }
      }
    }
  }
  graphics::legend("topright", legend = attributes, text.col = line.col, bty = "n", ncol = 2, text.font = legend.font,
                   cex = legend.cex)
  if (save.as != "") grDevices::dev.off()
}

#' Convert TCATA data
#'
#' Converts TCATA data from a set of onset-offset times to an indicator vector (\code{0}s and \code{1}s). Also works for TDS data.
#' @name convert.tcata
#' @aliases convert.tcata
#' @param X matrix with onset (start) times in first column and offset (stop) times in second column
#' @param times time slices for output indicator vector
#' @param decimal.places decimal places used in \code{times}; used for naming of the indices of \code{out.vec}
#' @return out.vec indictor vector(\code{0}s and \code{1}s)
#' @export
#' @encoding UTF-8
#' @examples
#' X <- rbind(c(3.18, 6.83), c(8.46, 11.09), c(18.61, 21.80))
#' times <- seq(0, 25, by = 0.01)
#' Xnew <- convert.tcata(X, times)
#' Xnew
convert.tcata <- function(X, times, decimal.places = 2){
  # matrix format:
  no.error <- TRUE
  if(ncol(X) != 2) no.error <- FALSE
  if(nrow(X) == 0) no.error <- FALSE
  array.base = NA
  if(times[1] == 0) array.base = 0
  if(times[1] == 1) array.base = 1
  if(is.na(array.base)) no.error <- FALSE
  adj <- ifelse(array.base == 0, 1, 0)
  #if(array.base == 0){
  #  adj <- 1
  #} else {
  #  adj <- 0
  #}
  times <- unique(round(times * (10^decimal.places), 0) / (10^decimal.places))
  out.vec <- rep(0, times = length(times))
  in.vec <- X * (10^decimal.places)

  if (no.error) {
    for( r in 1:nrow(in.vec) ){
      out.vec[ (in.vec[r, 1]+adj):(in.vec[r, 2]+adj) ] <- 1
    }
  } else {
    print( "convert.tcata() gives an error." )
  }
  return( out.vec )
}

#' Convert Temporal Category data
#'
#' Converts Temporal Category data from a set of onset-offset times and ratings to an vector of ratings.
#' @name convert.tcategory
#' @aliases convert.tcategory
#' @param X matrix with onset (start) times in first column and offset (stop) times in second column
#' @param in.scores vector of category values corresponding to rows of \code{X}
#' @param times time slices for output vector
#' @param decimal.places decimal places used in \code{times}; used for naming of the indices of \code{out.vec}
#' @return out.vec indictor vector(\code{0}s and \code{1}s)
#' @export
#' @encoding UTF-8
#' @examples
#' X <- rbind(c(3.18, 6.83), c(8.46, 11.09), c(18.61, 21.80))
#' in.scores <- c(7, 6, 5)
#' times <- seq(0, 25, by = 0.01)
#' Xnew <- convert.tcategory(X, in.scores, times)
#' Xnew
convert.tcategory <- function(X, in.scores, times, decimal.places = 2){
  # matrix format:
  no.error <- TRUE
  if(ncol(X) != 2) no.error <- FALSE
  if(nrow(X) == 0) no.error <- FALSE
  array.base = NA
  if(times[1] == 0) array.base = 0
  if(times[1] == 1) array.base = 1
  if(is.na(array.base)) no.error <- FALSE
  adj <- ifelse(array.base == 0, 1, 0)
  times <- unique(round(times * (10^decimal.places), 0) / (10^decimal.places))
  out.vec <- rep(0, times = length(times))
  in.vec <- X * (10^decimal.places)
  if(length(in.scores) != nrow(X)) no.error <- FALSE

  if(no.error) {
    for( r in 1:nrow(in.vec) ){
      out.vec[ (in.vec[r, 1]+adj):(in.vec[r, 2]+adj) ] <- in.scores[r]
    }
  } else {
    print( "convert.tcategory() gives an error." )
  }
  return( out.vec )
}

#' Get a pretty palette of colours
#'
#' Create a vector of n pretty colours.
#' @name pretty_palette
#' @aliases pretty_palette
#' @encoding UTF-8
#' @param n number of colours in the palette
#' @return cv A character vector, \code{cv}, of colours that look pretty.
#' @export
#' @examples
#' pretty_palette(8)
pretty_palette <- function(n){
  more.col = c()
  if(n > 12) {
    requireNamespace("grDevices", quietly = TRUE)
    more.col = grDevices::rainbow(n - 8, start = 0, end = 0.995)
  }
  return(c("red", "lightblue", "forestgreen", "purple", "hotpink", "chocolate4",
      "orange", "steelblue", "grey", "green",  "khaki", "maroon",
      more.col)[1:n])
}

#' Convenience function for getting a pretty palette and highlight colours
#'
#' Make a vector of n pretty colours, and n matching highlight colours.
#' @name make.palettes
#' @aliases make.palettes
#' @encoding UTF-8
#' @param n number of colours for each palette
#' @return pal A character vector, \code{cv}, of colours that look pretty.
#' @return pal.light A character vector, \code{cv}, of matching highlight colours that look pretty.
#' @export
#' @examples
#' make.palettes(8)
make.palettes <- function(n){
  pal <- pal.light <- pretty_palette(n)
  lighten.palette <- function(x){
    requireNamespace("grDevices", quietly = TRUE)
    x <- adjust.brightness(grDevices::col2rgb(x), 0)
    x <- adjust.brightness(grDevices::col2rgb(x), 2500)
    return(x)
  }
  pal.light <- sapply(pal, lighten.palette)
  names(pal.light) <- NULL
  return(list(pal = pal, pal.light = pal.light))
}

#' Plot trajectories based on Temporal Check-All-That-Apply (TCATA) data
#'
#' Plot trajectories following PCA on multiblock TCATA proportions, or same for Temporal Dominance of Sensations (TDS) proportions.
#' @aliases plot_pca.trajectories
#' @usage plot_pca.trajectories(in.pca = in.pca, products.times = matrix(NA),
#' attributes = c(), type = "smooth", span = 0.75, biplot = "distance",
#' flip = c(FALSE, FALSE), dims = c(1, 2),
#' att.offset.x = c(), att.offset.y = c(), att.cex = 1, inflate.factor = NA,
#' xlab = "_auto_", ylab = "_auto_", xlim = NULL, ylim = NULL,
#' attributes.col = "red", attributes.pch = 17,
#' lwd = 1, traj.lab.loc = 0, traj.col = c(grDevices::grey(1/2)), traj.points = NA,
#' traj.col.seg = NA, traj.cex = 1, traj.lab = c(), traj.lab.cex = 1,
#' arrow.loc = NA, arrow.length = 0.1, arrow.col = NA, arrow.lwd = NA,
#' contrails = list(), main = "", save.format = "eps", save.as = "")
#' @param in.pca Any \code{list} object with components \code{sdev}, \code{rotation}, and \code{x}. Most often it is a \code{prcomp} object obtained from PCA on a matrix of proportions (or, if there is no missing data, on counts) with Product*Times in rows and Attributes in columns.
#' @param products.times a 2-column matrix, with an ascending sort order on products (column 1) and a secondary ascending sort on times (column 2), corresponding to the rows of the matrix submitted to prcomp to obtain \code{"in.pca"}.
#' @param attributes a vector of attribute labels, corresponding to the attributes of the matrix submitted to prcomp to obtain \code{"in.pca"}.
#' @param type Determines how trajectories are drawn. Possible values are \code{"smooth"} (default) or \code{"raw"}.
#' @param span A tuning parameter used if smoothing trajectories using the \code{loess} function.
#' @param biplot Controls the type of biplot displayed. Possible values are \code{"distance"} (plots trajectories based on scores, and attributes based on eigenvectors multiplied by the \code{"inflation.factor"}), or \code{"correlation"} (plots trajectories based on scores divided by the sqrt of their respective eigenvalues, and attributes based on eigenvectors multiplied by the sqrt of their respective eigenvalues).
#' @param flip a vector of two logical values. Value indicates whether to mirror the coordinates in the x and y dimensions respectively. Default is \code{c(FALSE, FALSE)}.
#' @param dims a vector of two integers, specifying the principal componts to display. Defaults is \code{"c(1, 2)"}, i.e. PC1 vs. PC2.
#' @param att.offset.x A vector of numeric values corresponding to the labels in \code{"attributes"}. Used to adjust the horizontal position of attribute labels to make the plot more readable.
#' @param att.offset.y A vector of numeric values corresponding to the labels in \code{"attributes"}. Used to adjust the vertical position of attribute labels to make the plot more readable.
#' @param att.cex Attribute text size.
#' @param inflate.factor Scalar controlling the position of attribute labels. If \code{"NA"} (default), then this scalar is set to the largest absolute score divided by the largest absolute eigenvector based on the dimensions used. Use \code{"1"} for no inflation. Applies only when \code{biplot = "distance"}.
#' @param xlab Label for x axis.
#' @param ylab Label for y axis.
#' @param xlim Permits control of the x limit. Limits can be specified using a vector of 2 (ascending) numbers. If a single number is provided then values are selected such that the limits are 20\% beyond the smallest and largest x coordinates, respectively. If unspecified then control over x axis limits is given to the plot function in R.
#' @param ylim Permits control of the x limit using the same logic as is used for \code{"xlim"}.
#' @param attributes.col Color used to display attribute labels (see \code{"attributes"}).
#' @param attributes.pch Symbol for attribute coordinates.
#' @param lwd Trajectory line width.
#' @param traj.col A vector of colors for trajectories. If not specified then all trajectories are shown in grey.
#' @param traj.points Specifies the position of markers along smoothed trajectories, and used to indicate the progression of time.
#' @param traj.col.seg A vector of colors for segments along trajectories. If \code{NA} (default) then no segments along the trajectories appear in a color other than those specified by \code{"traj.col"}. This parameter applies to smoothed trajectories only.
#' @param traj.cex Used with \code{"traj.points"} for smoothed trajectories. Controls the size of symbol displayed.
#' @param traj.lab A vector of character labels that identify the trajectories. If unspecified, then products are identified by ascending natural numbers.
#' @param traj.lab.loc Indicates where along the trajectory the trajectory label will be positioned. \code{"1"} indicates the start of the trajectory. The value \code{"0"} (default) is a special convention indicating the end of the trajectory.
#' @param traj.lab.cex Text size of \code{traj.lab}.
#' @param arrow.loc Trajectory arrows locations for direction marker(s).
#' @param arrow.length Trajectory arrows length. See \code{length} parameter in \code{\link[graphics]{arrows}}.
#' @param arrow.col Trajectory arrows color. See \code{col} parameter in \code{\link[graphics]{arrows}}.
#' @param arrow.lwd Trajectory arrows line width. See \code{lwd} parameter in \code{\link[graphics]{arrows}}.
#' @param contrails list of data.frame objects with columns x, y, count, col; x and y are coordinates, count is the number of values at the coordinate, and col is the rbg colour.
#' @param main plot title; see \code{\link[graphics]{plot}}.
#' @param save.format If indicated, this will be the file type for the save image. Defaults to \code{"eps"} (eps format). Other possible values are \code{""} (not saved) or \code{"png"} (png format).
#' @param save.as The filename. Must be provided if the file will be saved.
#' @export
#' @seealso \code{\link[stats]{prcomp}}, \code{\link[graphics]{par}}
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @references Castura, J.C., Baker, A.K., & Ross, C.F. (2016). Using contrails and animated sequences to visualize uncertainty in dynamic sensory profiles obtained from temporal check-all-that-apply (TCATA) data. \emph{Food Quality and Preference}, 54, 90-100. \doi{10.1016/j.foodqual.2016.06.011}
#' @examples
#' # example using 'syrah' data set
#' syrah.pca <- prcomp(syrah[1:248, -c(1:4)], scale. = FALSE)
#' plot_pca.trajectories(syrah.pca, products.times = syrah[1:124, c(1, 4)],
#'                       attributes = colnames(syrah)[-c(1:4)], type = "raw")
#'
#' # now with smoothing; may need to play with the span parameter to get appropriate smoothing
#' plot_pca.trajectories(syrah.pca, products.times = syrah[1:124, c(1, 4)],
#'                       attributes = colnames(syrah)[-c(1:4)], type = "smooth", span = 0.3)
#'
#' # plots at each time point (trajectories join 2 points so start at timepoint 2, i.e., 11 s)
#' x <- 11:14 # for brevity show only the first 4 timeslices
#' # x <- 11:41 # uncomment this line to to run a longer demo
#' pca.list <- list()
#' for(i in seq_along(x)){
#'   pca.list[[x[i]-10]] <- syrah.pca
#'   pca.list[[x[i]-10]]$x <- pca.list[[x[i]-10]]$x[1:((x[i]-9)*6), ]
#'   plot_pca.trajectories(pca.list[[x[i]-10]], products.times = syrah[1:((x[i]-9)*6), c(1, 4)],
#'                         attributes = colnames(syrah)[-c(1:4)], type = "raw", inflate.factor = 1.5)
#'   Sys.sleep(3/4)
#'   # save plot if saving stills for a video; see Castura, Baker, & Ross (2016, Video 1)
#' }
#'
plot_pca.trajectories <- function( in.pca = in.pca,
                                   products.times = matrix(NA),
                                   attributes = c(),
                                   type = "smooth",
                                   span = 0.75,
                                   biplot = "distance",
                                   flip = c(FALSE, FALSE),
                                   dims = c(1, 2),
                                   att.offset.x = c(),
                                   att.offset.y = c(),
                                   att.cex = 1,
                                   inflate.factor = NA,
                                   xlab = "_auto_",
                                   ylab = "_auto_",
                                   xlim = NULL,
                                   ylim = NULL,
                                   attributes.col = "red",
                                   attributes.pch = 17,
                                   lwd = 1,
                                   traj.lab.loc = 0,
                                   traj.col = c(grDevices::grey(1/2)),
                                   traj.points = NA,
                                   traj.col.seg = NA,
                                   traj.cex = 1,
                                   traj.lab = c(),
                                   traj.lab.cex = 1,
                                   arrow.loc = NA,
                                   arrow.length = 0.1,
                                   arrow.col = NA,
                                   arrow.lwd = NA,
                                   contrails = list(),
                                   main = "",
                                   save.format = "eps",
                                   save.as = "" ){
  requireNamespace("grDevices", quietly = TRUE)
  requireNamespace("graphics", quietly = TRUE)
  # attribute loadings
  loadings <- in.pca$rotation[,dims] * in.pca$sdev[dims]
  if( length(attributes) > 0 ) {
   y.names <- attributes
  } else {
   y.names <- paste0("var",c(1:nrow(loadings)))
  }
  # rescale for biplot
  if(biplot == "distance"){
    scores <- in.pca$x[, dims]  # scores (products)
    maxscore <- max(abs(scores))
    eigenvectors <- in.pca$rotation[, dims]
    maxeigenvectors <- max(abs(eigenvectors))
    if(is.na(inflate.factor)){
      inflate.factor <- maxscore / maxeigenvectors
    }
    y <- eigenvectors * inflate.factor # new loadings (attributes)
  }
  if(biplot == "correlation"){
    scores <- in.pca$x[, dims]  # scores (products)
    scores[, 1] <- scores[, 1] / in.pca$sdev[dims[1]]
    scores[, 2] <- scores[, 2] / in.pca$sdev[dims[2]]
    y <- eigenvectors <- in.pca$rotation[,dims]
    y[, 1] <- y[, 1] * in.pca$sdev[dims[1]]
    y[, 2] <- y[, 2] * in.pca$sdev[dims[2]]
  }
  if(length(xlim) == 1 & is.numeric(xlim)){
    # treat xlim as a multiplier
    xlim <- range(c(scores[, 1], y[, 1])) * xlim
  }
  if(length(ylim) == 1 & is.numeric(ylim)){
    # treat xlim as a multiplier
    ylim <- range(c(scores[, 2], y[, 2])) * ylim
  }
  if(flip[1] == TRUE) {
    y[, 1] <- y[, 1]*(-1)
    scores[, 1] <- scores[, 1]*(-1)
  }
  if(flip[2] == TRUE) {
    y[, 2] <- y[, 2]*(-1)
    scores[, 2] <- scores[, 2]*(-1)
  }
  if(dim(products.times)[1]==1 & dim(products.times)[2]==1){
    print("Products & Time matrix missing. Defaulting to one product over time.")
    #** products.times <- cbind(1, 1:length(in.pca$x[, 1]))
    products.times <- data.frame(products = 1, times = 1:length(in.pca$x[, 1]))
  }
  products.labels <- unique(products.times[, 1])
  times.labels <- unique(products.times[, 2])
  if(length(att.offset.x) == 0) att.offset.x <- rep(0, nrow(in.pca$rotation))
  if(length(att.offset.y) == 0) att.offset.y <- rep(0, nrow(in.pca$rotation))

  if(length(att.offset.x) != nrow(in.pca$rotation) | length(att.offset.y) != nrow(in.pca$rotation))
    return(print(paste( "att.offset.x and att.offset.y must equal the number of attributes (",as.character( nrow( in.pca$rotation ) ), ")" ) ))
  # draw empty plot
  var.exp <- 100 * in.pca$sdev[dims]^2 / sum(in.pca$sdev^2)
  if (xlab == "_auto_")
    xlab = paste0("Dimension ", dims[1], " (", format(var.exp[1], nsmall = 2,  digits = 2), "%)")
  if (ylab == "_auto_")
    ylab = paste0("Dimension ", dims[2], " (", format(var.exp[2], nsmall = 2,  digits = 2), "%)")
  for (gg in 1:length(save.format)){
    if(save.format[gg] == "eps" & save.as[gg] != "") grDevices::postscript(save.as[gg])
    graphics::plot(c(y[, 1], scores[, 1]), c(y[, 2], scores[, 2]), type="n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, axes = FALSE, main = main)
    graphics::box()
    graphics::abline(h = 0, v = 0, lty = 3)
    graphics::axis(1)
    graphics::axis(2)
    DO_CONTRAILS <- FALSE
    if(length(contrails) > 0){
      if(length(contrails) == length(products.labels)){
        DO_CONTRAILS <- TRUE
        for(i in length(contrails)){
          # ensure the format of the contrails list is good - should have 1 data frame per product, with these columns: 'x', 'y', 'count'
          if(!all(c("x", "y", "count", "col") %in% colnames(contrails[[i]]))) DO_CONTRAILS <- FALSE
        }
      }
    }
    # add product trajectories
    for(p in 1:length(products.labels) ){
      this.product.scores <- scores[ products.times[, 1] == products.labels[p], ]
      if(DO_CONTRAILS){
        graphics::points(x = contrails[[p]]$x, y = contrails[[p]]$y, col = contrails[[p]]$col, pch = 16, cex = .5)
      }

      if(type=="raw" || length(this.product.scores[, 1]) < 10){
        for(this.time in 1:nrow(this.product.scores)){
          if(this.time > 1){
            last.time = this.time - 1
            graphics::lines(this.product.scores[last.time:this.time, 1], this.product.scores[last.time:this.time, 2],
                   lwd = lwd, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]))
          }
        }
        graphics::points(x = this.product.scores[traj.points, 1],
               y = this.product.scores[traj.points, 2],
               pch = 20, cex = traj.cex, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]) ) # plot trajectories
        if(!is.na(arrow.loc[1])){
          # Add arrow(s)
          if(is.na(arrow.col[1])) arrow.col <- traj.col
          if(is.na(arrow.lwd[1])) arrow.lwd <- lwd
          for(arr in seq_along(arrow.loc)){
            this.arrow.xy <- this.product.scores[arrow.loc[arr], 1:2]
            this.arrow.ang <- atan2(data.matrix(mean(this.product.scores[max(arrow.loc[arr], 1):min(arrow.loc[arr] + 25,
                                                                                                    length(this.product.scores)), 2])) - data.matrix(this.arrow.xy[2]),
                                    data.matrix(mean(this.product.scores[max(arrow.loc[arr], 1):min(arrow.loc[arr] + 25,
                                                                                                    length(this.product.scores)), 1])) - data.matrix(this.arrow.xy[1])) * 180 / pi
            # get the 'next point' inelegantly
            if(this.arrow.ang > 0) this.arrow.next.xy <- this.arrow.xy + c(1, 1)*arrow.length/10
            if(this.arrow.ang > 90) this.arrow.next.xy <- this.arrow.xy + c(-1, 1)*arrow.length/10
            if(this.arrow.ang < 0) this.arrow.next.xy <- this.arrow.xy + c(1, -1)*arrow.length/10
            if(this.arrow.ang < -90) this.arrow.next.xy <- this.arrow.xy + c(-1, -1)*arrow.length/10
            graphics::arrows(x0 = data.matrix(this.arrow.xy[1]), y0 = data.matrix(this.arrow.xy[2]),
                             x1 = data.matrix(this.arrow.next.xy[1]), y1 = data.matrix(this.arrow.next.xy[2]),
                             length = arrow.length, code = 2,
                             col = arrow.col[ifelse(length(arrow.col) > 1, p, 1)],
                             lwd = arrow.lwd[ifelse(length(arrow.lwd) > 1, p, 1)])
          }
        }
        # set location of the terminal product labels (end unless validly specified otherwise)
        if (traj.lab.loc == 0 | traj.lab.loc > length(this.product.scores[, 1])) {
          traj.lab.loc = length(this.product.scores[, 1])
        }
        graphics::points(x=this.product.scores[, 1][traj.lab.loc],
               y=this.product.scores[, 2][traj.lab.loc],
               pch = 15, cex = 2.75, lwd = lwd, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]) )
        graphics::text(x = this.product.scores[, 1][traj.lab.loc],
             y = this.product.scores[, 2][traj.lab.loc],
             labels = ifelse(length(traj.lab) == 0, as.character(p), as.character(traj.lab[p])), cex = traj.lab.cex, col="white")
      }
      if(type=="smooth" & length(this.product.scores[, 1]) >= 10){
        requireNamespace("stats", quietly = TRUE)
        x.loess <- stats::loess(this.product.scores[, 1] ~ I(1:nrow(this.product.scores)), span = span) #assume departures from smoothness are due to error
        y.loess <- stats::loess(this.product.scores[, 2] ~ I(1:nrow(this.product.scores)), span = span)
        graphics::points(x = stats::fitted(x.loess)[traj.points], y = stats::fitted(y.loess)[traj.points],
               pch = 20, cex=traj.cex, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]) ) # plot trajectories

        #if (class(traj.col.seg) == "logical" & is.na(traj.col.seg)[1]) {

        if(all(is.logical(traj.col.seg)) & is.na(traj.col.seg)[1]){
          graphics::lines(stats::fitted(x.loess), stats::fitted(y.loess), lwd = lwd, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]))
        } else {
          for(j in 2:length(stats::fitted(x.loess))){
            graphics::lines(stats::fitted(x.loess)[c(j - 1, j)],
                            stats::fitted(y.loess)[c(j - 1, j)],
                            lwd = lwd, col = traj.col.seg[p, j - 1])
          }
        }
        if(!is.na(arrow.loc[1])){
          # Add arrow(s)
          if(is.na(arrow.col[1])) arrow.col <- traj.col
          if(is.na(arrow.lwd[1])) arrow.lwd <- lwd
          for(arr in seq_along(arrow.loc)){
            this.arrow.xy <- c(stats::fitted(x.loess)[arrow.loc[arr]], stats::fitted(y.loess)[arrow.loc[arr]])
            this.arrow.ang <- atan2(mean(stats::fitted(y.loess)[max(arrow.loc[arr], 1):min(arrow.loc[arr] + 25, length(stats::fitted(y.loess)))]) - this.arrow.xy[2],
                                    mean(stats::fitted(x.loess)[max(arrow.loc[arr], 1):min(arrow.loc[arr] + 25, length(stats::fitted(x.loess)))]) - this.arrow.xy[1]) * 180 / pi
            # get the 'next point' inelegantly
            if(this.arrow.ang > 0) this.arrow.next.xy <- this.arrow.xy + c(1, 1)*arrow.length/100
            if(this.arrow.ang > 90) this.arrow.next.xy <- this.arrow.xy + c(-1, 1)*arrow.length/100
            if(this.arrow.ang < 0) this.arrow.next.xy <- this.arrow.xy + c(1, -1)*arrow.length/100
            if(this.arrow.ang < -90) this.arrow.next.xy <- this.arrow.xy + c(-1, -1)*arrow.length/100
            graphics::arrows(this.arrow.xy[1], this.arrow.xy[2], this.arrow.next.xy[1], this.arrow.next.xy[2],
                   length = arrow.length, code = 2,
                   col = arrow.col[ifelse(length(arrow.col) > 1, p, 1)],
                   lwd = arrow.lwd[ifelse(length(arrow.lwd) > 1, p, 1)])
          }
        }
        # set location of the terminal product labels (end unless validly specified otherwise)
        if (traj.lab.loc == 0 | traj.lab.loc > length(stats::fitted(x.loess))) {
          traj.lab.loc = length(stats::fitted(x.loess))
        }
        graphics::points(x = stats::fitted(x.loess)[traj.lab.loc], y = stats::fitted(y.loess)[traj.lab.loc],
               pch = 15, cex = 2.75, lwd = lwd, col = ifelse(length(traj.col) == 1, traj.col[1], traj.col[p]))
        graphics::text(x = stats::fitted(x.loess)[traj.lab.loc], y = stats::fitted(y.loess)[traj.lab.loc],
             labels=ifelse(length(traj.lab) == 0, as.character(p), as.character(traj.lab[p])), cex = traj.lab.cex, col="white")
      }
    }
    graphics::points(y[, 1], y[, 2], cex = 1, col = attributes.col, pch = attributes.pch)
    yoff1 <- .5 * graphics::strwidth(y.names, cex = 0.75) + .5 * graphics::strwidth("o", cex = .75)
    yoff2 <- .5 * graphics::strheight(y.names, cex = 0.75) + .5 * graphics::strheight("o", cex = .75)
    graphics::text(y[, 1] + yoff1 + att.offset.x, y[, 2] + yoff2 + att.offset.y, y.names, cex = att.cex, xpd = TRUE, col = attributes.col)

    if(save.format[gg] == "png" & save.as[gg] != "") grDevices::savePlot(save.as[gg], type = "png")
    if(save.format[gg] == "eps" & save.as[gg] != "") grDevices::dev.off()
  }
}

#' Calculate city block distance between two matrices
#'
#' Calculates the city block distance between two matrices.
#' @name dist.city.block
#' @aliases dist.city.block
#' @usage dist.city.block(x, y)
#' @param x first matrix
#' @param y second matrix
#' @return cbdist city block distance between \code{x} and \code{y}
#' @export
#' @encoding UTF-8
#' @examples
#'   x <- matrix(0, nrow = 5, ncol = 7)
#'   y <- matrix(1, nrow = 5, ncol = 7)
#'   dist.city.block(x, y)
#'
#'   y <- matrix(c(rep(0, 15), rep(1, 20)), nrow = 5, ncol = 7)
#'   dist.city.block(x, y)
dist.city.block <- function (x, y) {
  if(dim(x)[[1L]] == dim(y)[[1L]] & dim(x)[[2L]] == dim(y)[[2L]]){
    return(mean(abs(x - y), na.rm=TRUE))
  } else {
    print("Incompatible matrices in dist.city.block")
  }
}

#' Quantify TCATA assessor replication
#'
#' Quantify TCATA assessor replication using city block distance
#' @name similarity.tcata.replication
#' @aliases similarity.tcata.replication
#' @param this.assessor TCATA data (given as an indicator matrix) for assessor of interest
#' @param other.assessors TCATA data (given as an indicator matrix) for other assessors
#' @return replication.index city block distance between this assessor and other assessors
#' @details Similarity between one TCATA assessor and other assessors on the panel is quantified. The replication index can take on values between \code{0} and \code{1}, which indicate complete dissimilarity (disagreement) and complete similarity (agreement), respectively.
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @examples
#'   # Toy TCATA data for three assessors: a1, a2, a3
#'   a1 <- rbind(rep(0, 7),
#'               rep(0, 7),
#'               c(0, 0, 0, 1, 1, 1, 1),
#'               c(0, 0, 0, 1, 1, 1, 1),
#'               c(0, 0, 0, 1, 1, 1, 0))
#'   a2 <- rbind(c(0, 0, 0, 1, 1, 1, 0),
#'               rep(0, 7),
#'               c(0, 1, 1, 1, 1, 1, 0),
#'               rep(1, 7),
#'               c(0, 0, 0, 1, 1, 1, 1))
#'   a3 <- rbind(rep(0, 7),
#'               rep(0, 7),
#'               rep(1, 7),
#'               rep(1, 7),
#'               rep(1, 7))
#'
#'   # Quantify similarity of assessor a1 to the other assessors
#'   similarity.tcata.replication(a1, rbind(a2, a3))
similarity.tcata.replication <- function(this.assessor, other.assessors){
  n.blocks <- nrow(other.assessors) / nrow(this.assessor)
  if(n.blocks %% 1 != 0) print("Incompatible matrices in similarity.tcata.replication")
  if(dim(this.assessor)[[1L]] > dim(other.assessors)[[1L]] &
     dim(this.assessor)[[2L]] != dim(other.assessors)[[2L]]){
    print("Incompatible matrices in similarity.tcata.replication")
  }
  this.asssessor.block <- this.assessor
  if(n.blocks > 1){
    for(i in 2:n.blocks){
      this.asssessor.block <- rbind(this.asssessor.block, this.assessor)
    }
  }
  return(1-dist.city.block(this.asssessor.block, other.assessors))
}

#' Quantify TCATA assessor repeatability
#'
#' Quantify TCATA assessor repeatability using city block distance
#' @name similarity.tcata.repeatability
#' @aliases similarity.tcata.repeatability
#' @param X list of matrices, where each matrix is a TCATA data (given as an indicator matrix) for assessor of interest for one rep
#' @return repeatability.index average city block distance between matrices from replicated evaluations
#' @details Similarity between repeated evaluations given by a TCATA assessor is quantified. The repeatability index can take on values between \code{0} and \code{1}, which indicate complete dissimilarity (non-repeatability) and complete similarity (repeatability), respectively.
#' @export
#' @encoding UTF-8
#' @references Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
#' @examples
#'   # Toy data from one TCATA assessor on a product over three sessions: rep1, rep2, rep3
#'   rep1 <- rbind(rep(0, 7),
#'                 rep(0, 7),
#'                 c(0, 0, 0, 1, 1, 1, 1),
#'                 c(0, 0, 0, 1, 1, 1, 1),
#'                 c(0, 0, 0, 1, 1, 1, 0))
#'   rep2 <- rbind(c(0, 0, 0, 1, 1, 1, 0),
#'                 rep(0, 7),
#'                 c(0, 1, 1, 1, 1, 1, 0),
#'                 rep(1, 7),
#'                 c(0, 0, 0, 1, 1, 1, 1))
#'   rep3 <- rbind(rep(0, 7),
#'                 rep(0, 7),
#'                 rep(1, 7),
#'                 rep(1, 7),
#'                 rep(1, 7))
#'   rep.data <- list(rep1, rep2, rep3)
#'
#'   # Quantify similarity of assessor a1 to the other assessors
#'   similarity.tcata.repeatability(rep.data)
similarity.tcata.repeatability <- function(X){
  if(is.list(X)){
    for(ll in 2:length(X)){
      if(dim(X[[ll - 1]])[[1L]] != dim(X[[ll]])[[1L]] ||
         dim(X[[ll - 1]])[[2L]] != dim(X[[ll]])[[2L]]){
        print("Incompatible matrices in similarity.tcata.repeatability")
        return(NA)
      }
    }
  } else {
    return(1) # The agreement of something with itself is 1
  }
  requireNamespace("utils", quietly = TRUE)
  combn.reps <- utils::combn(length(X), 2)
  dis.repeat <- 0
  for(cc in 1:ncol(combn.reps)){
    dis.repeat <- dis.repeat + dist.city.block(X[[combn.reps[1, cc]]],
                                               X[[combn.reps[2, cc]]])
  }
  return(1-dis.repeat/ncol(combn.reps))
}

#' Count attribute selections
#'
#' Count the number of times that the attribute was selected (or optionally: deselected) in a single TCATA or TDS evaluation.
#' @name count.selections
#' @aliases count.selections
#' @param x vector of binary data (with possible values \code{0} or \code{1})
#' @param deselections set to \code{TRUE} if purpose is to count the number of deselections
#' @return count of selections (or deselections if \code{deselections = TRUE})
#' @details Count the number of times that the attribute was selected (or, optionally, deselected) in a single TCATA or TDS evaluation.
#' @export
#' @encoding UTF-8
#' @examples
#' data(bars)
#' paste0(bars[1, -c(1:4)], collapse = "")
#' # this attribute was checked 3 times and unchecked 2 times
#' count.selections(bars[1, -c(1:4)])
#' count.selections(bars[1, -c(1:4)], deselections = TRUE)
count.selections <- function(x, deselections = FALSE){
  search_pattern <- "01"
  append_pattern <- "0"
  if(deselections == TRUE){
    search_pattern <- "10"
    append_pattern <- "1"
  }
  x <- paste0("0", paste0(x, collapse = ""), append_pattern) # will not interfere with finding the search pattern
  return(length(unlist(strsplit(x, search_pattern, fixed = TRUE))) - 1)
}

#' Get bootstrap confidence bands for attribute selections
#'
#' Get bootstrap confidence bands for TCATA attribute citation rates or TDS attribute dominance rates.
#' @name bootstrap.band
#' @aliases bootstrap.band
#' @param X data frame of indicator data (with possible values \code{0} or \code{1})
#' @param boot number of virtual panels
#' @param alpha alpha level for bootstrap confidence bands
#' @param return.bias indicates whether to return bias associated with bootstrap mean value
#' @return \code{lcl} lower \code{100(alpha/2)\%} bootstrap confidence limit
#' @return \code{ccl} upper \code{100(1 - alpha/2)\%} bootstrap confidence limit
#' @return \code{bias} provided if \code{output.bias = TRUE}
#' @details Get bootstrap confidence bands for TCATA attribute citation rates or TDS attribute dominance rates.
#' @export
#' @encoding UTF-8
#' @examples
#' x <- ojtcata[ojtcata$samp == 1 & ojtcata$attribute == "Sweetness",  -c(1:4)]
#' x.boot.ci <- bootstrap.band(x, boot = 99) # 99 is only for illustrative purposes
#' x.boot.ci
bootstrap.band <- function(X, boot = 999, alpha = 0.05, return.bias = FALSE){
  tmp <- matrix(0, ncol = ncol(X), nrow = boot + 1)
  for (i in 1:boot){
    tmp[i, ] <- apply(X[sample(1:nrow(X), nrow(X), replace = TRUE), ], 2, mean)
  }
  tmp[boot + 1, ] <- apply(X, 2, mean) # last row is real data
  out <- data.frame(bias = apply(tmp, 2, mean) - apply(X, 2, mean), lcl = 0, ucl = 0)
  out[, 2:3] <- t(apply(tmp, 2, stats::quantile, probs = c(alpha/2, 1 - alpha/2))) # naive bands
  if(return.bias){
    return(list(lcl = out$lcl, ucl = out$ucl, bias = out$bias))
  } else {
    return(list(lcl = out$lcl, ucl = out$ucl))
  }
}



#' Draw h-cross, range box, and box to enclose h-box
#'
#' Draw h-cross, range box, and box to enclose h-cross, described
#' by Castura, Rutledge, Ross & Næs (2022).
#' @name draw.hcross
#' @aliases draw.hcross
#' @usage draw.hcross(rangebox = NULL, hcross = NULL,
#' rbox.col = "black", rbox.lty = "dotted", rbox.lwd = 4.5,
#' hbox.col = "lightgrey", hbox.lty = "solid", hbox.lwd = 4.5,
#' hcross.col = "black",hcross.lty = "solid",
#' hcross.signif.lwd = 7, hcross.nsd.lwd = 3.5)
#' @param rangebox matrix where columns 1 and 2 are x and y dimensions and
#' rows 1 and 2 are the minimum and maximum values
#' @param hcross matrix where columns 1 and 2 are x and y dimensions and
#' rows 1 and 2 are the half-width of the confidence interval, which is often
#' 95\% thus approximately 2x the standard error
#' @param rbox.col  line color for the range box (default: \code{"black"})
#' @param rbox.lty line type for the range box (default: \code{"dotted"})
#' @param rbox.lwd  line width for the range box (default: \code{4.5})
#' @param hbox.col  line color for the box enclosing the h-cross
#' (default: \code{"lightgrey"})
#' @param hbox.lty line type for the box enclosing the h-cross
#' (default: \code{"dotted"})
#' @param hbox.lwd  line width for the box enclosing the h-cross
#' (default: \code{4.5})
#' @param hcross.col line color for the h-cross (default: \code{"solid"})
#' @param hcross.lty line type for the h-cross (default: \code{"solid"})
#' @param hcross.signif.lwd line width for the h-cross where there is a
#' significant difference (default: \code{7})
#' @param hcross.nsd.lwd line width for the h-cross where there is a
#' significant difference (default: \code{3.5})
#' @details Draw h-cross, range box, and box to enclose h-box.
#' @references Castura, J.C., Rutledge, D.N., Ross, C.F., & Næs, T. (2022). Discriminability and uncertainty in principal component analysis (PCA) of temporal check-all-that-apply (TCATA) data. \emph{Food Quality and Preference}, 96, 104370. \doi{10.1016/j.foodqual.2021.104370}
#' @encoding UTF-8
draw.hcross <- function(rangebox = NULL, hcross = NULL,
                        rbox.col = "black", rbox.lty = "dotted", rbox.lwd = 4.5,
                        hbox.col = "lightgrey", hbox.lty = "solid", hbox.lwd = 4.5,
                        hcross.col = "black", hcross.lty = "solid",
                        hcross.signif.lwd = 7, hcross.nsd.lwd = 3.5){
  if(!is.null(rangebox[1])){
    # range box
    graphics::lines(rep(rangebox[1,1],2), rangebox[,2],
          lty = rbox.lty, lwd = rbox.lwd, col = rbox.col)
    graphics::lines(rep(rangebox[2,1],2), rangebox[,2],
          lty = rbox.lty, lwd = rbox.lwd, col = rbox.col)
    graphics::lines(rangebox[,1], rep(rangebox[1,2],2),
          lty = rbox.lty, lwd = rbox.lwd, col = rbox.col)
    graphics::lines(rangebox[,1], rep(rangebox[2,2],2),
          lty = rbox.lty, lwd = rbox.lwd, col = rbox.col)
  }
  if(!is.null(hcross[1])){
    # box for h-cross
    graphics::lines(rep(hcross[1,1],2), hcross[,2],
          lty = hbox.lty, lwd = hbox.lwd, col = hbox.col)
    graphics::lines(rep(hcross[2,1],2), hcross[,2],
          lty = hbox.lty, lwd = hbox.lwd, col = hbox.col)
    graphics::lines(hcross[,1], rep(hcross[1,2],2),
          lty = hbox.lty, lwd = hbox.lwd, col = hbox.col)
    graphics::lines(hcross[,1], rep(hcross[2,2],2),
          lty = hbox.lty, lwd = hbox.lwd, col = hbox.col)
  }
  if(!is.null(rangebox[1]) & !is.null(hcross[1])){
    # h-cross
    centerpoint <- apply(rangebox, 2, mean)
    graphics::lines(rep(centerpoint[1], 2), hcross[,2], lty = hcross.lty, col = hcross.col,
          lwd = ifelse(max(rangebox[,1])-min(rangebox[,1])-
                         max(hcross[,1])+min(hcross[,1]) > 0,
                       hcross.signif.lwd, hcross.nsd.lwd))
    graphics::lines(hcross[,1], rep(centerpoint[2], 2), lty = hcross.lty, col = hcross.col,
          lwd = ifelse(max(rangebox[,2])-min(rangebox[,2])-
                         max(hcross[,2])+min(hcross[,2]) > 0,
                       hcross.signif.lwd, hcross.nsd.lwd))
  }
  invisible()
}


