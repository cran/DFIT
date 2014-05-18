#e###############################################################################
# # powerExample.R
# # R Versions: 2.15.0
# #
# # Author(s): Victor H. Cervantes
# #
# # Description: Shows how to use the package DFIT to calculate power
# #
# # Outputs: Example runs for power calculation in the DFIT approach for NCDIF
# #
# # File history:
# #   20120403: Creation
################################################################################

################################################################################
# # Load packages
################################################################################
library(DFIT)
#source("../R/IPR.R")
#source("../R/IRTSE.R")
#source("../R/MantelHaenszel.R")
#source("../R/RajuAreas.R")
#source("../R/NCDIF.R")
#library(ggplot2)

################################################################################
# # Load functions
################################################################################

################################################################################
# # Definition of input and output paths
################################################################################
#inPath  <- "../input/"
#outPath <- "../output/"
#srcPath <- "../src/"
#logPath <- "../log/"

################################################################################
# # Other global variables
################################################################################
set.seed(45627)

nReplicates <- 500
nFocal      <- 500
nReference  <- 2500
kRatio      <- nReference / nFocal

################################################################################
# # Uniform power example
################################################################################

# # Dichotomous item parameters
dichotomousItemParameters <- list(focal = cbind(rep(1, 51), seq(-0.5, 0.5, length = 51)),
                                  reference = cbind(rep(1, 51), rep(0, 51)))

# # Two parameter logistic model normal metric
twoPlUniformParameters <- dichotomousItemParameters

twoPlUniformAse <- list()
twoPlUniformAse[["focal"]] <- AseIrt(itemParameters = twoPlUniformParameters[["focal"]],
                                     sampleSize = nFocal, irtModel = "2pl", logistic = FALSE)
twoPlUniformAse[["reference"]] <- AseIrt(itemParameters = twoPlUniformParameters[["reference"]],
                                         sampleSize = nReference, irtModel = "2pl", logistic = FALSE)



twoPlUniformNcdif    <- Ncdif(twoPlUniformParameters, logistic = FALSE)
twoPlUniformIpr      <- Ipr(itemParameters = twoPlUniformParameters, itemCovariances = twoPlUniformAse, nReplicates = nReplicates)
twoPlUniformNcdifIpr <- IprNcdif(itemParameterList = twoPlUniformIpr, irtModel = "2pl", logistic = FALSE, subdivisions = 1000)

# Calculate cutoff point
nullFocal                <- which(twoPlUniformParameters[['focal']][, 2] == twoPlUniformParameters[['reference']][, 2])
cutoffPointEachSZUniform <- CutoffIpr(quantiles = 0.95, iprStatistics = matrix(twoPlUniformNcdifIpr[nullFocal, ],
                                                                               nrow = length(nullFocal)))

# # Calculate cutoff with focal sample size
twoPlUniformParametersCutOff                <- twoPlUniformParameters
twoPlUniformParametersCutOff[['focal']]     <- matrix(twoPlUniformParametersCutOff[['focal']][nullFocal, ], nrow = length(nullFocal))
twoPlUniformParametersCutOff[["reference"]] <- matrix(twoPlUniformParametersCutOff[["reference"]][nullFocal, ], nrow = length(nullFocal))

twoPlUniformAseCutOff                <- twoPlUniformAse
twoPlUniformAseCutOff[['focal']]     <- twoPlUniformAseCutOff[['focal']][nullFocal]
twoPlUniformAseCutOff[['reference']] <- lapply(twoPlUniformAseCutOff[['reference']][nullFocal], '*', kRatio)

twoPlUniformIprCutOff <- Ipr(itemParameters = twoPlUniformParametersCutOff, itemCovariances = twoPlUniformAseCutOff, nReplicates = nReplicates)

twoPlUniformNcdifIprCutOff <- IprNcdif(itemParameterList = twoPlUniformIprCutOff, irtModel = "2pl", logistic = FALSE, subdivisions = 1000)
cutoffPointUniform         <- CutoffIpr(quantiles = 0.95, iprStatistics = matrix(twoPlUniformNcdifIprCutOff,
                                                                                 nrow = length(nullFocal)))



# Calculate power
powerUniform <- data.frame(b = twoPlUniformParameters[['focal']][, 2],
                           NCDIF = twoPlUniformNcdif *
                                   sign(twoPlUniformParameters[['focal']][, 2] - 0),
                           Power = rowMeans(twoPlUniformNcdifIpr > cutoffPointUniform$quantiles))

powerPlotUniform <- ggplot(powerUniform, aes(x = b, y = Power))
powerPlotUniform <- powerPlotUniform + geom_line()

powerPlotUniform
# ggsave(paste(outPath, 'uniformPower.pdf', sep = ''), width = 8, height = 5)

powerPlotUniform <- ggplot(powerUniform, aes(x = NCDIF, y = Power))
powerPlotUniform <- powerPlotUniform + geom_line()

powerPlotUniform
#ggsave(paste(outPath, 'uniformPowerNCDIF.pdf', sep = ''), width = 8, height = 5)


powerUniformEach <- data.frame(Type = "Cutoff with both groups",
                               b = twoPlUniformParameters[['focal']][, 2],
                               NCDIF = twoPlUniformNcdif *
                                       sign(twoPlUniformParameters[['focal']][, 2] - 0),
                               Power = rowMeans(twoPlUniformNcdifIpr > cutoffPointEachSZUniform$quantiles))

powerUniform <- data.frame(Type = "Cutoff with only focal group",
                           powerUniform)

powerUniform <- rbind(powerUniform, powerUniformEach)

powerPlotUniform <- ggplot(powerUniform, aes(x = NCDIF, y = Power))
powerPlotUniform <- powerPlotUniform + geom_line(aes(colour = Type, linetype = Type), size = 1.2)
powerPlotUniform <- powerPlotUniform + xlim(-0.015, 0.015)

powerPlotUniform
#ggsave(paste(outPath, 'uniformPowerNCDIFEach.pdf', sep = ''), width = 10, height = 5)




################################################################################
# # Non uniform power example
################################################################################

# # Dichotomous item parameters
dichotomousItemParameters <- list(focal = cbind(unique(c(seq(0.5, 1, length = 26), 1, seq(1, 1 / 0.5, length = 26))), rep(0, 51)),
                                  reference = cbind(rep(1, 51), rep(0, 51)))


# # Two parameter logistic model normal metric
twoPlNonUniformParameters <- dichotomousItemParameters

twoPlNonUniformAse                <- list()
twoPlNonUniformAse[["focal"]]     <- AseIrt(itemParameters = twoPlNonUniformParameters[["focal"]],
                                            sampleSize = nFocal, irtModel = "2pl", logistic = FALSE)
twoPlNonUniformAse[["reference"]] <- AseIrt(itemParameters = twoPlNonUniformParameters[["reference"]],
                                            sampleSize = nReference, irtModel = "2pl", logistic = FALSE)



twoPlNonUniformNcdif    <- Ncdif(twoPlNonUniformParameters, logistic = FALSE)
twoPlNonUniformIpr      <- Ipr(itemParameters = twoPlNonUniformParameters, itemCovariances = twoPlNonUniformAse, nReplicates = nReplicates)
twoPlNonUniformNcdifIpr <- IprNcdif(itemParameterList = twoPlNonUniformIpr, irtModel = "2pl", logistic = FALSE, subdivisions = 1000)

# Calculate cutoff point
nullFocal                   <- which(twoPlNonUniformParameters[['focal']][, 1] == twoPlNonUniformParameters[['reference']][, 1])
cutoffPointEachSZNonUniform <- CutoffIpr(quantiles = 0.95, iprStatistics = matrix(twoPlNonUniformNcdifIpr[nullFocal, ],
                                                                                  nrow = length(nullFocal)))

# # Calculate cutoff with focal sample size
twoPlNonUniformParametersCutOff                <- twoPlNonUniformParameters
twoPlNonUniformParametersCutOff[["reference"]] <- matrix(twoPlNonUniformParametersCutOff[["reference"]][nullFocal, ], nrow = length(nullFocal))
twoPlNonUniformParametersCutOff[['focal']]     <- matrix(twoPlNonUniformParametersCutOff[['focal']][nullFocal, ], nrow = length(nullFocal))

twoPlNonUniformAseCutOff                <- twoPlNonUniformAse
twoPlNonUniformAseCutOff[["reference"]] <- lapply(twoPlNonUniformAseCutOff[["reference"]][nullFocal],  '*',  kRatio)
twoPlNonUniformAseCutOff[['focal']]     <- twoPlNonUniformAseCutOff[['focal']][nullFocal]

twoPlNonUniformIprCutOff <- Ipr(itemParameters = twoPlNonUniformParametersCutOff, itemCovariances = twoPlNonUniformAseCutOff, nReplicates = nReplicates)

twoPlNonUniformNcdifIprCutOff <- IprNcdif(itemParameterList = twoPlNonUniformIprCutOff, irtModel = "2pl", logistic = FALSE, subdivisions = 1000)
cutoffPointNonUniform         <- CutoffIpr(quantiles = 0.95, iprStatistics = matrix(twoPlNonUniformNcdifIprCutOff,
                                                                                    nrow = length(nullFocal)))


# Calculate power
powerNonUniform <- data.frame(a = twoPlNonUniformParameters[['focal']][, 1],
                              NCDIF = twoPlNonUniformNcdif *
                                      sign(twoPlNonUniformParameters[['focal']][, 1] - 1),
                              Power = rowMeans(twoPlNonUniformNcdifIpr > cutoffPointNonUniform$quantiles))

powerPlotNonUniform <- ggplot(powerNonUniform, aes(x = a, y = Power))
powerPlotNonUniform <- powerPlotNonUniform + geom_line()

powerPlotNonUniform
#ggsave(paste(outPath, 'nonuniformPower.pdf', sep = ''), width = 8, height = 5)

powerPlotNonUniform <- ggplot(powerNonUniform, aes(x = NCDIF, y = Power))
powerPlotNonUniform <- powerPlotNonUniform + geom_line()

powerPlotNonUniform
#ggsave(paste(outPath, 'nonuniformPowerNCDIF.pdf', sep = ''), width = 8, height = 5)


powerNonUniformEach <- data.frame(Type = "Cutoff with both groups",
                                  a = twoPlNonUniformParameters[['focal']][, 1],
                                  NCDIF = twoPlNonUniformNcdif *
                                          sign(twoPlNonUniformParameters[['focal']][, 1] - 1),
                                  Power = rowMeans(twoPlNonUniformNcdifIpr > cutoffPointEachSZNonUniform$quantiles))

powerNonUniform <- data.frame(Type = "Cutoff with only focal group",
                              powerNonUniform)

powerNonUniform <- rbind(powerNonUniform, powerNonUniformEach)

powerPlotNonUniform <- ggplot(powerNonUniform, aes(x = NCDIF, y = Power))
powerPlotNonUniform <- powerPlotNonUniform + geom_line(aes(colour = Type, linetype = Type), size = 1.2)
powerPlotNonUniform <- powerPlotNonUniform + xlim(-0.015, 0.015)

powerPlotNonUniform
#ggsave(paste(outPath, 'nonuniformPowerNCDIFEach.pdf', sep = ''), width = 10, height = 5)



# # Calculate bias
biasUniform <- data.frame(NCDIF = twoPlUniformNcdif *
                                  sign(twoPlUniformParameters[['focal']][, 2] - 0),
                          Bias = rowMeans(twoPlUniformNcdifIpr) - twoPlUniformNcdif)

biasUniform[, 'RelBias'] <- abs(biasUniform[, 'Bias'] / biasUniform[, 'NCDIF'])
biasUniform[, 'Type']    <- 'Uniform DIF'

biasNonUniform <- data.frame(NCDIF = twoPlNonUniformNcdif *
                                     sign(twoPlNonUniformParameters[['focal']][, 1] - 1),
                             Bias = rowMeans(twoPlNonUniformNcdifIpr) - twoPlNonUniformNcdif)

biasNonUniform[, 'RelBias'] <- abs(biasNonUniform[, 'Bias'] / biasUniform[, 'NCDIF'])
biasNonUniform[, 'Type']    <- 'Nonuniform DIF'

biasUniform <- rbind(biasUniform, biasNonUniform)

biasPlotUniform <- ggplot(biasUniform, aes(x = NCDIF, y = Bias))
biasPlotUniform <- biasPlotUniform + geom_line(aes(colour = Type)) + ylab('Bias') + xlim(-0.015, 0.015) + ylim(0, 0.015)

biasPlotUniform
#ggsave(paste(outPath, 'uniformBias.pdf', sep = ''), width = 8, height = 5)

relbiasPlotUniform <- ggplot(biasUniform, aes(x = NCDIF, y = RelBias))
relbiasPlotUniform <- relbiasPlotUniform + geom_line(aes(colour = Type)) + ylab('Relative bias') + xlim(-0.015, 0.015)

relbiasPlotUniform
#ggsave(paste(outPath, 'uniformRelativeBias.pdf', sep = ''), width = 8, height = 5)

