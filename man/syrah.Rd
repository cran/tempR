\alias{syrah}
\name{TCATA data set: Syrah wines}
\docType{data}
\title{TCATA data set: Syrah wines}
\description{TCATA citation proportions for three wine treatments evaluated using a two-sip evaluation protocol.
}
\format{
A data frame with 1026 rows (3 treatments * 2 sips * 171 time slices) and 13 columns:
    \itemize{
      \item [, 1] WineSip (chr) Code for wine and sip
      \item [, 2] Wine (chr) Code for wine (H=high, L=low, A=adjusted)
      \item [, 3] Sip (int) Sip number
      \item [, 4] Time (int) Time, in seconds
      \item [, 5] Astringency (num) citation proportions
      \item [, 6] Bitter (num) citation proportions
      \item [, 7] Dark Fruit (num) citation proportions
      \item [, 8] Earthy (num) citation proportions
      \item [, 9] Green (num) citation proportions
      \item [,10] Heat (num) citation proportions
      \item [,11] Red Fruit (num) citation proportions
      \item [,12] Spices (num) citation proportions
      \item [,13] Sour (num) citation proportions
    }
}

\encoding{UTF-8}
\references{Baker, A. K., Castura, J. C., & Ross, C. F. (2016). Temporal check-all-that-apply characterization of Syrah wine finish. \emph{Journal of Food Science}, 81, S1521-S1529. \url{http://dx.doi.org/10.1111/1750-3841.13328}.}
\examples{
head(syrah, 3) # review first 3 rows of 'syrah' data set
}
\keyword{ datasets }
\keyword{ TCATA }
