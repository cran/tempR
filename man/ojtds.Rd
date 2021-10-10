\alias{ojtds}
\name{TDS data set: orange juice}
\docType{data}
\title{TDS data set: orange juice}
\description{Raw results from 20-s TDS evaluations of six orange juice samples by 50 consumers.}
\format{
A data frame with 1800 rows (50 consumers * 6 samples * 6 attributes) and 25 columns (4 headers + 21 time slices)
\itemize{
  \item [, 1] cons (int) consumer id
  \item [, 2] samp (chr) sample id
  \item [, 3] samp_pos (int) position of sample in serving order
  \item [, 4] attribute (chr) sensory attribute
  \item [, 5:25]  time_\code{99}s (int) value is \code{1} if attribute is selected at time slice; otherwise value is \code{0}
  }
}
\encoding{UTF-8}
\references{Ares, G., Jaeger, S. R., Antúnez, L., Vidal, L, Giménez, A., Coste, B., Picallo, A., & Castura, J.C. (2016). Comparison of TCATA and TDS for dynamic sensory characterization of food products. \emph{Food Research International}, 78, 148-158. \doi{10.1016/j.foodres.2015.10.023}}
\examples{
head(ojtds) # review first 6 rows of 'ojtds' data set
}
\keyword{ datasets }
\keyword{ TDS }
