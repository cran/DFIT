% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPR.R
\name{Bound3PlIpr}
\alias{Bound3PlIpr}
\title{Takes item parameters from Ipr and forces guessing to lie between 0 and 1}
\usage{
Bound3PlIpr(itemParameterList)
}
\arguments{
\item{itemParameterList}{A list where each element is a list containing "focal" and "reference" item Parameters from a 3PL model. Item parameters are assumed to be on the same scale. Item parameters for each group should be a matrix with nrow equal to the number of items}
}
\value{
itemParameterList     A list where each element is a list containing "focal" and "reference" item Parameters where guessing parameters outside the [0, 1] interval are changed by 0 and 1.
}
\description{
Makes all simulated guessing values from a 3PL model that are outside the [0, 1] interval to be 0 or 1.
}
\examples{
# # Not run
# #
# # data(dichotomousItemParameters)
# # threePlParameters <- dichotomousItemParameters
# # isNot3Pl          <- ((dichotomousItemParameters[['focal']][, 3] == 0) |
# #                      (dichotomousItemParameters[['reference']][, 3] == 0))
# #
# # threePlParameters[['focal']]          <- threePlParameters[['focal']][!isNot3Pl, ]
# # threePlParameters[['reference']]      <- threePlParameters[['reference']][!isNot3Pl, ]
# # threePlParameters[['focal']][, 3]     <- threePlParameters[['focal']][, 3] + 0.1
# # threePlParameters[['reference']][, 3] <- threePlParameters[['reference']][, 3] + 0.1
# # threePlParameters[['focal']][, 2]     <- threePlParameters[['focal']][, 2] + 1.5
# # threePlParameters[['reference']][, 2] <- threePlParameters[['reference']][, 2] + 1.5
# # threePlParameters[['focal']]          <- threePlParameters[['focal']][-c(12, 16, 28), ]
# # threePlParameters[['reference']]      <- threePlParameters[['reference']][-c(12, 16, 28), ]
# #
# # threePlAse <- list()
# # threePlAse[["focal"]]     <- AseIrt(itemParameters = threePlParameters[["focal"]],
# #                                     logistic = TRUE,
# #                                     sampleSize = 10000,
# #                                     irtModel = "3pl")
# # threePlAse[["reference"]] <- AseIrt(itemParameters = threePlParameters[["reference"]],
# #                                     logistic = TRUE,
# #                                     sampleSize = 15000,
# #                                     irtModel = "3pl")
# #
# # set.seed(41568)
# # threePlIpr <- Ipr(itemParameters = threePlParameters, itemCovariances = threePlAse,
# #                   nReplicates = 100)
# # threePlIpr <- Bound3PlIpr(threePlIpr)

}
\author{
Victor H. Cervantes <vhcervantesb at unal.edu.co>
}
