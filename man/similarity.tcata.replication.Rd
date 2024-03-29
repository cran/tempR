% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcata.R
\encoding{UTF-8}
\name{similarity.tcata.replication}
\alias{similarity.tcata.replication}
\title{Quantify TCATA assessor replication}
\usage{
similarity.tcata.replication(this.assessor, other.assessors)
}
\arguments{
\item{this.assessor}{TCATA data (given as an indicator matrix) for assessor of interest}

\item{other.assessors}{TCATA data (given as an indicator matrix) for other assessors}
}
\value{
replication.index city block distance between this assessor and other assessors
}
\description{
Quantify TCATA assessor replication using city block distance
}
\details{
Similarity between one TCATA assessor and other assessors on the panel is quantified. The replication index can take on values between \code{0} and \code{1}, which indicate complete dissimilarity (disagreement) and complete similarity (agreement), respectively.
}
\examples{
  # Toy TCATA data for three assessors: a1, a2, a3
  a1 <- rbind(rep(0, 7),
              rep(0, 7),
              c(0, 0, 0, 1, 1, 1, 1),
              c(0, 0, 0, 1, 1, 1, 1),
              c(0, 0, 0, 1, 1, 1, 0))
  a2 <- rbind(c(0, 0, 0, 1, 1, 1, 0),
              rep(0, 7),
              c(0, 1, 1, 1, 1, 1, 0),
              rep(1, 7),
              c(0, 0, 0, 1, 1, 1, 1))
  a3 <- rbind(rep(0, 7),
              rep(0, 7),
              rep(1, 7),
              rep(1, 7),
              rep(1, 7))

  # Quantify similarity of assessor a1 to the other assessors
  similarity.tcata.replication(a1, rbind(a2, a3))
}
\references{
Castura, J.C., Antúnez, L., Giménez, A., Ares, G. (2016). Temporal check-all-that-apply (TCATA): A novel temporal sensory method for characterizing products. \emph{Food Quality and Preference}, 47, 79-90. \doi{10.1016/j.foodqual.2015.06.017}
}
