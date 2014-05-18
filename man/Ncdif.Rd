\name{Ncdif}
\alias{Ncdif}
\title{Calculates NCDIF index for an item with given item parameters of focal and reference groups.}
\usage{
Ncdif(itemParameters, irtModel = "2pl", focalAbilities = NULL,
  focalDistribution = "norm", subdivisions = 5000, logistic = TRUE,
  focalDistrExtra = list(mean = 0, sd = 1))
}
\arguments{
  \item{itemParameters}{A list containing "focal" and
  "reference" item parameters. Item parameters are assumed
  to be on the same scale. Item parameters for each group
  should me a matrix with nrow equal to the number of
  items.}

  \item{irtModel}{A string stating the irtModel to be used.
  Should be one of "1pl", "2pl", "3pl", "grm" or "pcm".}

  \item{focalAbilities}{If NULL, NCDIF is calculated by
  numerical integration of focal distribution. If not NULL,
  it must be a numerical vector containing the abilities
  for the individuals in the focal group.}

  \item{focalDistribution}{A string stating the
  distribution name to be used for integrating. Only used
  if focalAbilities is NULL.}

  \item{focalDistrExtra}{Extra parameters for the focal
  group distribution function if needed.}

  \item{subdivisions}{A numeric value indicating the number
  of subdivisions for numerical integration. Only used if
  focalAbilities is NULL.}

  \item{logistic}{A logical value stating if the IRT model
  will use the logistic or the normal metric. Defaults to
  using the logistic metric by fixing the D constant to 1.
  If FALSE the constant is set to 1.702 so that the normal
  metric is used.}
}
\value{
ncdif Numeric vector with the NCDIF index value for each
item.
}
\description{
Calculates NCDIF index for an item with given item
parameters of focal and reference groups.
}
\examples{
data(dichotomousItemParameters)
threePlNcdif <- Ncdif(itemParameters = dichotomousItemParameters, irtModel = '3pl',
                      focalAbilities = NULL, focalDistribution = "norm",
                      subdivisions = 5000, logistic = TRUE)
}
\author{
Victor H. Cervantes <vcervantes at icfes.gov.co>
<vhcervantesb at unal.edu.co>
}
\references{
Raju, N. S., van der Linden, W. J., & Fleer, P. F. (1995).
An IRT-based internal measure of test bias with
applications for differential item functioning. Applied
Psychological Measurement, 19, 353--368.
}

