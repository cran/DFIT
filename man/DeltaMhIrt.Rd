\name{DeltaMhIrt}
\alias{DeltaMhIrt}
\title{Obtains the ETS Delta measure for Mantel-Haneszel DIF statistic effect size.}
\usage{
DeltaMhIrt(mh)
}
\arguments{
  \item{mh}{A numeric vector containing the MH statistic
  values}
}
\value{
delta A numeric vector containing the delta values
}
\description{
Obtains the ETS Delta measure for Mantel-Haneszel DIF
statistic effect size.
}
\examples{
data(dichotomousItemParameters)
threePlMh <- IrtMh(itemParameters = dichotomousItemParameters,  irtModel = "3pl",
                   focalDistribution = "norm",
                   referenceDistribution = "norm", focalDistrExtra = list(mean = 0, sd = 1),
                   referenceDistrExtra = list(mean = 0, sd = 1), groupRatio = 1,
                   logistic = FALSE)
delta3pl <- DeltaMhIrt(threePlMh)
}
\author{
Victor H. Cervantes <vcervantes at icfes.gov.co>
<vhcervantesb at unal.edu.co>
}
\references{
Holland, P.W., and Thayer, D.T. (1988). Differential Item
Performance and the Mantel-Haenszel Procedure. In H. Wainer
and H.I. Braun (Eds.), Test Validity. Hillsdale, NJ:
Erlbaum.
}

