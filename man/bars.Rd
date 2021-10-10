\alias{bars}
\name{TDS data set: snack bars}
\docType{data}
\title{TDS data set: snack bars}
\description{Raw TDS results from 24 assessors who evaluated four snack bars in triplicate.}
\format{
A data frame with 1440 rows (24 assessors * 3 sessions * 4 samples * 5 attributes) and 455 columns (4 header rows + 451 time slices)
\itemize{
  \item [,1] 	assessor (chr) assessor id
  \item [,2] 	session (chr) session id
  \item [,3] 	sample (chr) sample id
  \item [,4]  attribute (chr) sensory attribute
  \item [,5:455]  time_\code{99.9}s (chr) value is \code{1} if attribute is dominant at time slice; otherwise value is \code{0}
  }
}
\encoding{UTF-8}
\references{Findlay, C.J., Castura, J.C., & Valeriote, E. (2014). Temporal methods: A comparative study of four different techniques. In \emph{17th IUFoST Congress}. 17-21 August. Montréal, Québec, Canada.}
\examples{
head(bars, 2) # review first 2 rows of 'bars' data set
}
\keyword{ datasets }
\keyword{ TDS }
