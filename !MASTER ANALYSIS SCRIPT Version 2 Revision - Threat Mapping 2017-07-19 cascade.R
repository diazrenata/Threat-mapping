## MULTI-THREAT ASSESSMENT - THREAT MAPPING - MASTER ANALYSIS
## Adam B. Smith | Missouri Botanical Garden | adamDOTsmithATmobotDOTorg | January 2017
## Renata M. Diaz | Missouri Botanical Garden | diazDOTrenatamATgmailDOTcom | February 2017

# source('C:/ecology/Drive/Research/Multi-Threat Assessment/Analysis - Threat Mapping/!MASTER ANALYSIS SCRIPT Version 2 Revision - Threat Mapping 2017-01-18.r')



## Install necessary packages
# install.packages(c('sp', 'rgdal', 'raster', 'dismo','rgeos','geosphere','scales','RColorBrewer','TeachingDemos', 'doBy', 'phangorn', 'vegan', 'plyr', 'rgrass7','rstan','shinystan','RPushbullet'))


memory.limit(memory.limit() * 2^30)
rm(list=ls())
gc()

options(warnPartialMatchArgs=TRUE, warnPartialMatchAttr=TRUE, warnPartialMatchDollar=TRUE, stringsAsFactors=FALSE)

library(compiler)
enableJIT(1)
setCompilerOptions(suppressUndefined=TRUE)

library(sp)
library(rgdal)
library(raster)
options(java.parameters='-Xmx1g' )
library(dismo)
library(rgeos)
library(geosphere)
library(scales)
library(RColorBrewer)
library(TeachingDemos)
library(doBy)
library(phangorn)
library(vegan)
library(plyr)
library(shinystan)
library(rgrass7)
library(RPushbullet)
library(jsonlite)
library(slam)

cat('\n', date(), '\n'); flush.console()

################
### CONTENTS ###
################

### definitions ###

### BASIC ANALYSIS ###
### tally basic stats on threats ###
### barplots of proportions of species affected by each threat ###
### calculate percent threatened versus time of assessment by each threat ###
### histograms of threat frequency ###
### analyze relationship between number of threats and G-rank/number of counties occupied ###

### GEOGRAPHICAL ANALYSIS ###
### prepare GADM ###
### add number of species in each county by threat to GADM ###
### make maps of number of species with each threat type per county for species in just ONE county ###
### make maps of number of species and percentages of species with each threat type per county ###
### make maps of number of species and WITH and WITHOUT each threat type per county for species in just ONE county ###

### LITERATURE ANALYSIS ###
### literature analysis & plot ###

### THREAT SYNDROME ANALYSIS ###
### ordination ###
### hierarchical clustering analysis on THREATS ###
### hierarchical clustering analysis on SPECIES & THREATS ###
### calculate phi coefficient between threats ###

### THREAT MAPPING ANALYSIS ###
### tally number of counties in which each species occurs ###
### tally number of species found in each county
### map number of endemics per county ###
### plot histograms of date last assessed and date last updated ###


### THREAT MAPPING ANALYSIS ###
### Gather threats data ###
### add threat columns to natureServeCounties shapefile ###
### analysis with stan ###

#### Gather threats data ####
# load prefix function
# load natureServeCounties and add leading zeros to FIPS codes
# Climate change
# Climate stability
# Roads
# Human population density - 2010 census
# Venter (Human footprint) layers:
# Human footprint index - Venter
# Population density 1990 - Venter
# Population density 2010 - Venter
# Built environments - Venter
# Croplands - Venter
# Pasture - Venter
# Night lights - Venter
# Roads - Venter
# National land cover database
# Active mines ####
# MRDS
# Sriram et al layers:
# Sriram cars
# Sriram soil degradation
# Sriram proximity to roads
# Sriram pigs
# Sriram poultry
# Zabel crop suitability
# Human appropriation of NPP
# PADUS
# Invasive alien species - EDDMapS
# NLCD Percent Developed Imperviousness
# SPARROW Water quality (N and P pollution)
# Railroads
# Oil and gas wells
# Crude oil pipelines
# Shale plays
# Coal mines
#
#
# # load prefix function (redundant w/above)
#
# prefix <- function(x, len, pre='0') {
#
#   # x		value to add leading characters (will be coerced to character class)
#   # len	desired number of characters in x... will not be made shorter if x is longer than this
#   # pre	value to pre-pend to x, will be repeated until nchar(x)==len
#
#   x <- as.character(x)
#   size <- nchar(x)
#   if (nchar(x) < len) {
#     addTo <- paste(rep(pre, each=len - size), collapse='')
#     x <- paste0(addTo, x)
#   }
#   return(x)
# }
#
#
# # load natureServeCounties and add leading zeros to FIPS codes
# natureServeCounties <- shapefile('H:/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/ORIGINAL/NS_mv_CTY_bdrys_G12ESAtots_201403')
# natureServeCounties$FIPS_CODE_LONG <- apply(as.matrix(natureServeCounties$FIPS), 1, prefix, len = 5, pre = "0")



# ### climate change ###

# ## add 2014 and 2015 PRISM data to H-drive ##

#makeBilGTIFF <- function(bilFile) {
# if(substr(bilFile, as.numeric(nchar(bilFile)) - 3, as.numeric(nchar(bilFile))) == '.xml') return ()
# thisBil <- raster(paste0(workDir, '/', bilFile))
# fileTitle <- substr(bilFile, 0, as.numeric(nchar(bilFile)) - 4)
# writeRaster(thisBil, paste0(workDir, '/', fileTitle), format = 'GTiff', overwrite = TRUE, progress = 'text')
#}

#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/ppt/Recent years'
#bils <- as.matrix(list.files(workDir, full.names = FALSE, pattern = 'bil.bil'))
#apply(bils, 1, FUN = makeBilGTIFF)

#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmean/Recent years'
#bils <- as.matrix(list.files(workDir, full.names = FALSE, pattern = 'bil.bil'))
#apply(bils, 1, FUN = makeBilGTIFF)

#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmax/Recent years'
#bils <- as.matrix(list.files(workDir, full.names = FALSE, pattern = 'bil.bil'))
#apply(bils, 1, FUN = makeBilGTIFF)

#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmin/Recent years'
#bils <- as.matrix(list.files(workDir, full.names = FALSE, pattern = 'bil.bil'))
#apply(bils, 1, FUN = makeBilGTIFF)


#pptAnnualSum <- function(year){
#  months <- stack(list.files(workDir, full.names = TRUE, pattern = year))
#  sumRaster <- sum(months)
#  sumRaster <- setMinMax(sumRaster)
#  projection(sumRaster) <- projection(months)
#  writeRaster(sumRaster, paste0(workDir, '/ppt_annual_sum/ppt_sum_',year), format = 'GTiff', overwrite =  TRUE, progress = 'text')
#}

#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/ppt'
#pptAnnualSum('2014')
#pptAnnualSum('2015')


#tempAnnualMax <- function(year){
#  months <- stack(list.files(workDir, full.names = TRUE, pattern = year))
#  maxRaster <- max(months)
#  maxRaster <- setMinMax(maxRaster)
#  projection(maxRaster) <- projection(months)
#  writeRaster(maxRaster, paste0(workDir, '/tmax_annual_max/tmax_max_', year), format = 'GTiff', overwrite = TRUE, progress = 'text')
#}
#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmax'
#tempAnnualMax('2014')
#tempAnnualMax('2015')


#tempAnnualMean <- function(year){
#  months <- stack(list.files(workDir, full.names = TRUE, pattern = year))
#  meanRaster <- mean(months)
#  meanRaster <- setMinMax(meanRaster)
#  projection(meanRaster) <- projection(months)
#  writeRaster(meanRaster, paste0(workDir, '/tmean_annual_mean/tmean_mean_', year), format = 'GTiff', overwrite = TRUE, progress = 'text')
#}
#tempAnnualMean('2014')
#tempAnnualMean('2015')

#tempAnnualMin <- function(year){
#  months <- stack(list.files(workDir, full.names = TRUE, pattern = year))
#  minRaster <- min(months)
#  minRaster <- setMinMax(minRaster)
#  projection(minRaster) <- projection(months)
#  writeRaster(minRaster, paste0(workDir, '/tmin_annual_min/tmin_min_', year), format = 'GTiff', overwrite = TRUE, progress = 'text')
#}
#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmin'
#tempAnnualMin('2014')
#tempAnnualMin('2015')



# ## regressions for change in temp and precip ##

# ## precip ##
# run this once and then save .csv files to save time
# workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/ppt/ppt_annual_sum'
# usaStack <- stack(list.files(workDir, full.names = TRUE, pattern = "ppt_sum_"))
#
# # # need to make coordinates consistent
# # # PRISM data does not include Alaska or Hawaii or Puerto Rico
#
# natureServeCountiesClimate <- spTransform(natureServeCounties, CRS(projection(usaStack[[1]])))
# natureServeCountiesClimate <- natureServeCountiesClimate[ which(natureServeCountiesClimate$STATE_NAME != "Alaska"), ]
# natureServeCountiesClimate <- natureServeCountiesClimate[ which(natureServeCountiesClimate$STATE_NAME != "Hawaii"), ]
# natureServeCountiesClimate <- natureServeCountiesClimate[ which(natureServeCountiesClimate$STATE_NAME != "Puerto Rico"), ]
#
# getCountyAnnualMeans <- function(county_FIPS){
#   thisCounty <- natureServeCountiesClimate[ which(natureServeCountiesClimate$FIPS_CODE_LONG == county_FIPS), ]
#
#   years <- c(1895:2015)
#   annualMeans <- as.data.frame(years)
#
#   getAnnualMean <- function(year)  {
#     yearIndex <- as.numeric(match(year, annualMeans$years))
#     thisYear <- extract(usaStack[[yearIndex]], thisCounty, method = 'simple', na.rm = TRUE)
#     thisYearDf <- as.data.frame(thisYear)
#     thisYearMean <- mean(thisYearDf[ ,1], na.rm = TRUE)
#     return(thisYearMean)
#   }
#
#   annualMeans$countyAnnualMean <- apply(as.matrix(annualMeans$years), 1, getAnnualMean)
#
#   return(annualMeans)
# }
#
#
#
#
# getCountyRegressions <- function(x, climateVariable) {
#
#   if (climateVariable == "precip") thisCountyDf <- as.data.frame(precipAnnualMeans[x])
#   if (climateVariable == "tmean") thisCountyDf <- as.data.frame(tmeanAnnualMeans[x])
#   if (climateVariable == "precipProportional") {
#     thisCountyDf <- as.data.frame(precipAnnualMeans[x])
#     allYearsMean <- mean(thisCountyDf$countyAnnualMean)
#     thisCountyDf$countyAnnualMean <- thisCountyDf$countyAnnualMean / allYearsMean
#   }
#
#   if (is.na(thisCountyDf$countyAnnualMean[1])) return (c(NA, NA))
#   allYearsRegression <- lm(countyAnnualMean ~ years, data = thisCountyDf)
#   allYearsCoeff <- summary(allYearsRegression)$coefficients[2, 1]
#
#   pastThirty <- thisCountyDf[ which(thisCountyDf$years > 1984), ]
#   pastThirtyRegression <- lm(countyAnnualMean ~ years, data = pastThirty)
#   pastThirtyCoeff <- summary(pastThirtyRegression)$coefficients[2,1]
#   countyCoefficients <- c(allYearsCoeff, pastThirtyCoeff)
#
#
#   return(countyCoefficients)
# }
#
#
#
# precipAnnualMeans <- apply(as.matrix(natureServeCountiesClimate$FIPS_CODE_LONG), 1, getCountyAnnualMeans)
#
# precipRegressions <- apply(as.matrix(1:length(natureServeCountiesClimate$FIPS_CODE_LONG)), 1, getCountyRegressions, climateVariable = "precip")
#
# precipProportionalRegressions <- apply(as.matrix(1:length(natureServeCountiesClimate$FIPS_CODE_LONG)), 1, getCountyRegressions, climateVariable = 'precipProportional')
#
# natureServeCountiesClimate$precipAllYears <- precipRegressions[1, ]
# natureServeCountiesClimate$precipLastThirty <- precipRegressions[2, ]
#
# natureServeCountiesClimate$precipPropAllYears <- precipProportionalRegressions[1, ]
# natureServeCountiesClimate$precipPropLastThirty <- precipProportionalRegressions[2, ]
#
# rm(usaStack)
#
# ## temp ##
# workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmean/tmean_annual_mean'
# usaStack <- stack(list.files(workDir, full.names = TRUE, pattern = "tmean_mean_"))
#
# tmeanAnnualMeans <- apply(as.matrix(natureServeCountiesClimate$FIPS_CODE_LONG), 1, getCountyAnnualMeans)
# tmeanRegressions <- apply(as.matrix(1:length(natureServeCountiesClimate$FIPS_CODE_LONG)), 1, getCountyRegressions, climateVariable = "tmean")
#
# natureServeCountiesClimate$tmeanAllYears <- tmeanRegressions[1,]
# natureServeCountiesClimate$tmeanLastThirty <- tmeanRegressions[2,]
#
# rm(usaStack)
#
# ### Save this so you have the option of not running getCountyAnnualMeans.
# write.csv(natureServeCountiesClimate, 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/climate_change/regressions/natureServeCounties_withAllRegressions.csv')


# ##  Climate stability  ##
# # Characterize "past" and "present" climates for each county
# # Calculate bioclim variables for entire US for every year.
# yearlyBioclim <- function(year) {
#  tminDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmin'
#  tmaxDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/tmax'
#  precipDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/ppt'

#  tminMonthly <- stack(list.files(path = tminDir, pattern = paste('_', year, sep=""), full.names = TRUE))
#  tmaxMonthly <- stack(list.files(path = tmaxDir, pattern = paste('_', year, sep=""), full.names = TRUE))
#  precipMonthly <- stack(list.files(path = precipDir, pattern =paste('_', year, sep=""), full.names = TRUE))

#  bioClim <- biovars(precipMonthly, tminMonthly, tmaxMonthly)
#  saveDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual'

#  bioClimNames <- names(bioClim)

#  for (i in 1:19) {
#    writeRaster(bioClim[[i]], filename = paste(saveDir, "/", bioClimNames[i], "_", year, sep = ""), format = 'GTiff', overwrite = TRUE, progress = 'text')
#  }
#}

#years <- c(1895:2015)
#years <- as.character(years)
#years <- as.matrix(years)
#apply(years, 1, yearlyBioclim)

# # Get one dataframe of all variables for all pixels for all years
# # Columns are bio1, bio2, bio7, bio12, bio15, bio17.
# # Rows are cell-years from the entire-US rasters.

# getBioVarDf <- function(bioVariable) {
#  workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual'
#  setwd(workDir)
#  yearLayers <- stack(list.files(path = workDir, pattern = paste(bioVariable, '_', sep = ""), full.names=TRUE))

#  getYearDf <- function(yearIndex){
#    yearDf <- as.data.frame(yearLayers[[yearIndex]], na.rm = TRUE)
#    yearDf$pixel <- rownames(yearDf)
#    yearDf$year <- substr(colnames(yearDf)[1], nchar(colnames(yearDf)[1]) - 3, nchar(colnames(yearDf[1])))
#    colnames(yearDf) <- c(bioVariable, 'pixel','year')
#    return(yearDf)
#  }

#  for (j in 0:(floor(nlayers(yearLayers)/10))) {
#    startIndex <- j * 10 + 1
#    endIndex <- min((j + 1) * 10, nlayers(yearLayers))

#    variableDf <- getYearDf(startIndex)
#    if (endIndex > startIndex ) {
#      for (i in (startIndex + 1):endIndex){
#        variableDf <- rbind(variableDf, getYearDf(i))
#      }
#    }
#    workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/dfs'
#    setwd(workDir)
#    write.csv(variableDf, paste(bioVariable, "_start", startIndex, ".csv", sep = ""))
#    rm(variableDf)
#  }
#  rm(yearLayers)
#}


#getBioVarDf('bio1')
#getBioVarDf('bio2')
#getBioVarDf('bio7')
#getBioVarDf('bio12')
#getBioVarDf('bio15')
#getBioVarDf('bio17')

#bioVariable <- 'bio1'
#workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/dfs'
#segmentFiles <- list.files(path = workDir, pattern = paste(bioVariable, '_start', sep = ""), full.names=TRUE)
#bioVarDf <- read.csv(segmentFiles[1], header=TRUE)
#for (i in 2:length(segmentFiles)) {
# addingDf <- read.csv(segmentFiles[i], header = TRUE)
# bioVarDf <- rbind(bioVarDf, addingDf)
# rm(addingDf)
#}


#rm(bioVariable)
#rm(segmentFiles)
#rm(workDir)
#rm(i)
#bioVarDf <- bioVarDf[2:4]
#bioVarDf$bio2 <- NA
#bioVarDf$bio7 <- NA
#bioVarDf$bio12 <- NA
#bioVarDf$bio15 <- NA
#bioVarDf$bio17 <- NA

#addVariable <- function(bioVariable){
#  workDir <- 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/dfs'
#  segmentFiles <- list.files(path = workDir, pattern = paste(bioVariable, '_start', sep = ""), full.names=TRUE)
#  bioVarAdding <- read.csv(segmentFiles[1], header=TRUE)
#  for (i in 2:length(segmentFiles)) {
#    addingDf <- read.csv(segmentFiles[i], header = TRUE)
#    bioVarAdding <- rbind(bioVarAdding, addingDf)
#    rm(addingDf)
#  }
#  bioVarAdding <- bioVarAdding[,2]
#  return(bioVarAdding)
#}

#bioVar2 <- addVariable('bio2')
#bioVarDf$bio2 <- bioVar2
#rm(bioVar2)

#bioVar7 <- addVariable('bio7')
#bioVarDf$bio7 <- bioVar7
#rm(bioVar7)

#bioVar12 <- addVariable('bio12')
#bioVarDf$bio12 <- bioVar12
#rm(bioVar12)

#bioVar15 <- addVariable('bio15')
#bioVarDf$bio15 <- bioVar15
#rm(bioVar15)

#bioVar17 <- addVariable('bio17')
#bioVarDf$bio17 <- bioVar17
#rm(bioVar17)

#rm(addVariable)

#pcaCols <- c('bio1', 'bio2','bio7','bio12','bio15', 'bio17')
#bioVarDf <- bioVarDf[ , pcaCols]

#rm(pcaCols)
#gc()

#write.csv(bioVarDf, "bioVarDf.csv")

# # PCA on a random subset of all data

# ### loading only a portion of bioVarDf for troubleshooting.
# # rm(list=ls())
# bioVarDf <- read.csv('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/dfs/bioVarDf.csv')
# # bioVarDf <- read.csv('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/dfs/bioVarDf.csv', nrows = 100000)
#
# bioVarDf <- bioVarDf[,2:7]
# head(bioVarDf)
#
# set.seed(1)
# bioVarDfSub <- bioVarDf[ sample(nrow(bioVarDf), 10000000, replace = FALSE), ]
# bioVarDfSub <- bioVarDf[ sample(nrow(bioVarDf), 1000, replace = FALSE), ]

# head(bioVarDfSub)
# this is what head(bioVarDfSub) looks like if you load all rows, subsample 10000000 rows, and set.seed(1).
# bio1      bio2  bio7   bio12     bio15  bio17
# 29684381 13.217083 11.890833 37.95 1576.88  44.11326 237.33
# 55351384 12.655000 17.316667 34.70  576.39 102.29734  40.99
# 55582501  6.660000 11.471667 40.37 1230.78  68.07568 123.47
# 8035317  14.879167 11.780000 35.66 1591.76  53.04887 249.30
# 9434711  10.312083 14.912500 38.81  195.65  79.33290  17.74
# 22174615  1.425417  9.444167 32.91  851.09  61.05026  81.94

# rm(bioVarDf)

# pcaSub <- prcomp(bioVarDfSub, scale = TRUE)

# rm(bioVarDfSub)

# summary(pcaSub)
# Using the full sample:
# Importance of components:
#  PC1    PC2    PC3     PC4    PC5     PC6
#  Standard deviation     1.7386 1.1941 0.8430 0.62664 0.5618 0.36382
#  Proportion of Variance 0.5038 0.2376 0.1184 0.06545 0.0526 0.02206
#  Cumulative Proportion  0.5038 0.7415 0.8599 0.92533 0.9779 1.00000
#
#
# # Aggregate US raster and select every 4th year of data

# usaBio1 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio1_', full.names=TRUE))
# usaBio2 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio2_', full.names=TRUE))
# usaBio7 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio7_', full.names=TRUE))
# usaBio12 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio12_', full.names=TRUE))
# usaBio15 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio15_', full.names=TRUE))
# usaBio17 <- stack(list.files(path = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual', pattern = 'bio17_', full.names=TRUE))
#
# wholeUS <- stack(usaBio1, usaBio2, usaBio7, usaBio12, usaBio15, usaBio17)
# nlayers(wholeUS)

# The smallest county has 163 cells at highest resolution:
# natureServeCountiesClimate <- spTransform(natureServeCounties, CRS(proj4string(wholeUS)))
# smallestCty <- natureServeCountiesClimate[ which(natureServeCountiesClimate$Shape_Area == min(natureServeCountiesClimate$Shape_Area)), ]
# smallestCtyCells <- extract(wholeUS[[1]], smallestCty, cellnumbers = TRUE, df = TRUE)

# Rows and columns odd. Expand false.
# ncol(wholeUS[[1]])
# nrow(wholeUS[[1]])
#
# layers.to.keep <- 0:29
# layers.to.keep <- layers.to.keep * 4
# layers.to.keep <- layers.to.keep + 4
# layers.to.keep.start <- layers.to.keep
# for (i in 1:5) {
#   layers.to.keep <-c(layers.to.keep, layers.to.keep.start + (i * 121))
# }
#
# wholeUS.every4 <- wholeUS[[c(layers.to.keep)]]
#
# wholeUS.every4.ag <- aggregate(wholeUS.every4, fact = 2, fun = mean, expand = FALSE, na.rm = TRUE)
#
# # Create a df of all biovars for every 4th year for every aggregated cell
# cellsmax <- ncell(wholeUS.every4.ag[[1]])
# filled.cells <- extract(wholeUS.every4.ag[[1]], 1:cellsmax, df=TRUE)
#
# filled.cells <- na.omit(filled.cells)
# filled.cells <- filled.cells['ID']
# filled.cells <- as.vector(filled.cells[1])
# head(filled.cells)
#
# pullBioVars <- function(cellIndex) {
#   thisCell <- extract(wholeUS.every4.ag, filled.cells$ID[cellIndex], df = TRUE)
#   values <- as.matrix(thisCell, colnames = FALSE)
#   values <- values[2:181]
#
#   newMat <- matrix(data = values, nrow=30, ncol=6, byrow = FALSE, dimnames = NULL)
#   newMat <- as.data.frame(newMat)
#   colnames(newMat) <- c('bio1', 'bio2', 'bio7', 'bio12', 'bio15', 'bio17')
#   newMat$year <- layers.to.keep.start + 1894
#   newMat$cell <- filled.cells$ID[cellIndex]
#   write.csv(newMat, file = paste0('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/cells/every4ag_', filled.cells$ID[cellIndex], '.csv'), row.names = FALSE)
# }
#
# sapply(1:length(filled.cells$ID), pullBioVars)
#
#  # Use KDE to compare past and present climates
# ##  kde operates on every CELL across YEARS
#
# library(MASS)
#
# # cellIndex <- 1
#
# calculateKernelStability <- function(cellIndex){
#
#   thisCellAllYears <- read.csv(paste('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/cells/every4ag_', filled.cells$ID[cellIndex], '.csv', sep=""))
#   tail(thisCellAllYears)
#
#   thisCellAllYears <- thisCellAllYears[,1:6]
#   head(thisCellAllYears)
#
#   #### Get PCA component axes for two time periods ####
#   thisCellAllYears.pca <- predict(pcaSub, newdata = thisCellAllYears)
#   tail(thisCellAllYears.pca)
#
#
#   thisCellPastThirty <- thisCellAllYears[23:30,  ]
#   thisCellPastThirty.pca <- predict(pcaSub, newdata = thisCellPastThirty)
#   tail(thisCellPastThirty.pca)
#
#   #### Define past and present kernels ####
#   pastKernel <- kde2d(thisCellAllYears.pca[,1], thisCellAllYears.pca[,2], lims=c(range(thisCellAllYears.pca[,1]), range(thisCellAllYears.pca[,2])))
#   presentKernel <- kde2d(thisCellPastThirty.pca[,1], thisCellPastThirty.pca[,2], lims=c(range(thisCellAllYears.pca[,1]), range(thisCellAllYears.pca[,2])))
#
#   pastDf <- as.data.frame(pastKernel[3])
#   presentDf <- as.data.frame(presentKernel[3])
#
#   #### Standardize kernels so they add to 1 ####
#   pastDfStd <- pastDf/(sum(pastDf))
#   presentDfStd <- presentDf/(sum(presentDf))
#
#   #### Compute S ####
#   sumProbabilities <- 0
#   for (i in 1:25) {
#     for (j in 1:25) {
#       zPresent <- presentDfStd[i, j]
#       zPast <- pastDfStd[i, j]
#       d <- abs(zPast - zPresent)
#       sumProbabilities <- sumProbabilities + d
#     }
#   }
#
#   stability <- 1 - (sumProbabilities/2)
#
#   return(stability)
# }
#
# pixelStabilities <- lapply(1:length(filled.cells$ID), calculateKernelStability)
# pixelStabilities <- unlist(pixelStabilities)
# cellStabilities <- cbind(filled.cells, pixelStabilities)
#
# write.csv(cellStabilities, 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/climate_change/kde_stability/kdeStabilityIndex.csv')
# cellStabilities <- read.csv('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/climate_change/kde_stability/kdeStabilityIndex.csv')
#
# ## Make a raster of stabilities
#
# stability <- wholeUSevery4.ag[[1]]
# stability <- stability * 0 + 1
#
#
# stabilityfilled <- extract(stability, 1:ncell(stability), df=TRUE)
#
# stabilityfilled <- stabilityfilled[,1]
# stabilityfilled <- as.data.frame(stabilityfilled)
# colnames(stabilityfilled) <- 'ID'
# stabilityfilled <- join(stabilityfilled, cellStabilities, by = 'ID', type = 'left')
#
# stability <- setValues(stability, stabilityfilled$pixelStabilities, index = stabilityfilled$ID)
# names(stability) <- 'kde.Stability'
# plot(stability)
# stability
#
# writeRaster(stability, filename = 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/kde_stability', format = 'GTiff', overwrite = TRUE, progress = 'text')
#
# rm(cellStabilities)
# rm(stabilityfilled)
#
# #### Extract mean S per county ####
#
# kdeStability <- raster('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/kde_stability.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(kdeStability)))
# kdeStabilities <- natureServeCountiesLayer$FIPS_CODE_LONG
# kdeStabilities <- as.data.frame(kdeStabilities)
# colnames(kdeStabilities) <- 'FIPS_CODE_LONG'
# kdeStabilities$meanStability <- NA
# for (i in 1:nrow(kdeStabilities)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == kdeStabilities$FIPS_CODE_LONG[i]), ]
#   kdeStabilities$meanStability[i] <- extract(kdeStability, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# kdeStabilities <- as.matrix(kdeStabilities)
#
# write.csv(kdeStabilities, 'H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/kde_stabilities.csv')
#
#### Roads ####
#
#### Tally road lengths using GRASS ####
#
# states <- list.files(path ='H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/', pattern ='.shp', full.names=FALSE)
# clipstates <- function(statelong) {
#   statenew <- substr(statelong, (as.numeric(nchar(statelong)) - 5), as.numeric(nchar(statelong)) - 4)
#   return(statenew)
# }
# states<-sapply(states, clipstates)
#
# library(rgrass7)
#
# initGRASS('C:/Program Files/GRASS GIS 7.0.5', home=tempdir(), location = 'contus', mapset = 'PERMANENT', override = TRUE)
# get.GIS_LOCK()
# set.GIS_LOCK()
#
# getRoadLength <- function(state){
#   if (state == "ma"){
#     stateShape <- 'H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/road100k_l_ma.shp'
#     stateMapset <- 'roads_ma'
#     stateLayer <- 'road100k_l_ma'
#     stateLayerMapset <- paste(stateLayer, '@', stateMapset, sep = "")
#     stateDirectory <- 'H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/roads_ma'
#   } else {
#     stateShape <- paste('H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/street100k_l_', state, '.shp', sep='')
#     stateMapset <- paste('streets_', state, sep='')
#     stateLayer <- paste('street100k_l_', state, sep='')
#     stateLayerMapset <- paste(stateLayer, '@', stateMapset, sep="")
#     stateDirectory <- paste('H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/streets_', state, sep='')
#   }
#   execGRASS('g.mapset', mapset='PERMANENT', location='contus')
#   execGRASS('g.proj', flags='c', georef=stateShape)
#   execGRASS('g.mapset', flags='c', mapset='stateMapset', location='contus')
#   execGRASS('v.in.ogr', input=stateShape, layer=stateLayer, output=stateLayer)
#   execGRASS('v.db.addcolumn', map=stateLayer, columns='LENGTH double precision')
#   execGRASS('v.to.db', map=stateLayer, option='length', type='line', columns='LENGTH', units='meters')
#   execGRASS('db.out.ogr', input=stateLayer, output=stateDirectory, format='CSV')
# }
#
# for (i in 1:(as.numeric(length(states)))) {
#   print(states[i])
#   getRoadLength(states[i])
# }
#
# unset.GIS_LOCK()
# remove_GISRC()
# detach(package:rgrass7, unload = TRUE)

#
# #### Generate a df of road lengths per county ####
# ### Exclude MA, because MA only has roads (not streets)
# addState <- function(state){
#   if (state == 'ma'){
#     stateFile <- 'H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/roads_ma/roads_ma.csv'
#   } else {
#     stateFile <- paste('H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/streets_', state, '/streets_', state, '.csv', sep='')
#   }
#   stateDf <- read.csv(stateFile)
#   stateDf$STATEFP_LONG <- lapply(stateDf$STATEFP, prefix, len=2, pre='0')
#   stateDf$COUNTYFP_LONG <- lapply(stateDf$COUNTYFP, prefix, len=3, pre='0')
#   stateDf$FIPS_CODE_LONG <- paste(stateDf$STATEFP_LONG, stateDf$COUNTYFP_LONG, sep="")
#   colsKeep <- c('FIPS_CODE_LONG', 'cat', 'LENGTH')
#   stateDf <- stateDf[,colsKeep]
#   return(stateDf)
# }
#
# for (i in 1:(as.numeric(length(states)))) {
#   if(states[i] == 'ma') next
#   if (exists('allStatesDf')) {
#     allStatesDf <- rbind(allStatesDf, addState(states[i]))
#   } else {
#     allStatesDf <- addState(states[i])
#   }
# }
#
# countyLengths <- tapply(allStatesDf$LENGTH, allStatesDf$FIPS_CODE_LONG, sum)
# countyRoadLengths <- as.data.frame(countyLengths)
# countyRoadLengths$FIPS_CODE_LONG <- names(countyLengths)
# colnames(countyRoadLengths) <- c('LENGTH', 'FIPS_CODE_LONG')
# rm(countyLengths)
#
# # getwd()
# # setwd('H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/')
# # write.csv(countyRoadLengths, "county_road_lengths.csv")

#### Human population density - 2010 census ####
# path: H:/Global Change Program/GIS/Population/2010 Census Demographic Profile SF/DEC_10_DP_DPDP1_with_ann.csv
### Load and clean 2010 US census data; incorporate revisions
#
# workDir <- 'H:/Global Change Program/GIS/Population/2010 Census Demographic Profile SF'
# censusRaw <- read.csv(paste(workDir, '/DEC_10_DP_DPDP1_with_ann.csv', sep=""), skip = 1)
# censusRaw <- censusRaw[,1:4]
#
# getFIPS <- function(x) {
#   FIPS <- strsplit(as.character(x), "US")
#   FIPS <- unlist(FIPS)
#   FIPS <- FIPS[2]
# }
#
# censusRaw$FIPS_CODE_LONG <- sapply(censusRaw$Id, getFIPS)
# rm(getFIPS)
#
# cols <- colnames(censusRaw)
# cols[4] <- "Total_population"
# colnames(censusRaw) <- cols
# rm(cols)
#
# ## revisions
# # censusRevise <- censusRaw[which(censusRaw$Revisions == TRUE), ]
#
# revisionsInfo <- read.csv(paste(workDir, '/DEC_10_DP_DPDP1.txt', sep = ""), skip = 5, header = FALSE, sep= " ", blank.lines.skip = TRUE)
# revisionsInfo$rowIndex <- rownames(revisionsInfo)
# revisionIds <- as.data.frame(revisionsInfo[, c(1,11)])
#
# colnames(revisionIds) <- c('revisionId', 'rowIndex')
#
# findIds <- function(x){
#   isPar <- FALSE
#   if(substr(as.character(x), 1, 1) == '(') isPar <- TRUE
#   return(isPar)
# }
# revisionIds$keep <- sapply(revisionIds$revisionId, findIds)
#
# revisionIds <- revisionIds[which(revisionIds$keep == TRUE), 1:2]
#
# cleanIDs <- function(x){
#   id <- strsplit(as.character(x), ')')
#   id <- unlist(id)
#   id <- id[1]
#   id <- substr(id, 3, nchar(id))
#   return(id)
# }
#
# revisionIds$revisionId <- sapply(revisionIds$revisionId, cleanIDs)
# head(revisionIds)
#
# add1 <- function(x){
#   x <- as.numeric(x)
#   x <- x + 1
#   return(x)
# }
# revisionIds$valueIndex <- sapply(revisionIds$rowIndex, add1)
#
# valueIndices <- revisionIds$valueIndex
# valueIndices <- unlist(valueIndices)
#
# revisionValues <- as.data.frame(revisionsInfo[valueIndices, c(3, 11)])
# colnames(revisionValues) <- c('Revised_population', 'valueIndex')
#
# revisions <- cbind(revisionIds, revisionValues)
#
# revisions <- revisions[c('revisionId', 'Revised_population')]
#
# rm(revisionIds)
# rm(revisionsInfo)
# rm(valueIndices)
# rm(revisionValues)
# rm(add1)
# rm(cleanIDs)
# rm(findIds)
#
# censusRaw$revisionId <- NA
#
# cleanIDs.2 <- function(x){
#   if(grepl('r', as.character(x), fixed = FALSE) == FALSE) return(NA)
#   x <- strsplit(as.character(x), 'r')
#   x <- unlist(x)
#   x <- x[2]
#   x <- substr(x, 0, nchar(x) - 1)
#   return(x)
# }
#
# censusRaw$revisionId <- sapply(censusRaw$Total_population, cleanIDs.2)
# rm(cleanIDs.2)
#
# revisions$revisionId <- as.character(revisions$revisionId)
# censusRaw$revisionId <- as.character(censusRaw$revisionId)
#
# censusRevised <- join(censusRaw, revisions, by="revisionId", type = 'left')
#
# censusRevised$Revised_population <- gsub(",", censusRevised$Revised_population, replacement ="")
#
# joinPopulations <- function(x) {
#   if(is.na(censusRevised$Revised_population[x])) return(as.character(censusRevised$Total_population[x]))
#   return(as.character(censusRevised$Revised_population[x]))
# }
#
# censusRevised$Final_population <- sapply(c(1:length(censusRevised$FIPS_CODE_LONG)), joinPopulations)
# rm(joinPopulations)
#
# filepath <- 'H:/Global Change Program/GIS/Population/2010 Census Demographic Profile SF/censusClean.csv'
# write.csv(censusRevised, filepath)
# rm(revisions)
# rm(censusRaw)
# rm(filepath)

# ### Venter (Human footprint) layers ####
#
# extractVenterMean <- function(FIPS_CODE, venterLayer){
#   thisCounty <- natureServeCountiesVenter[ which(natureServeCountiesVenter$FIPS_CODE_LONG == FIPS_CODE), ]
#   thisCountyMean <- extract(venterLayer, thisCounty, fun = mean, na.rm=TRUE)
#   return(thisCountyMean)
# }
# #
# # extractVenterSum <- function(FIPS_CODE, venterLayer){
# #   thisCounty <- natureServeCountiesVenter[ which(natureServeCountiesVenter$FIPS_CODE_LONG == FIPS_CODE), ]
# #  thisCountySum <- extract(venterLayer, thisCounty, fun = sum, na.rm=TRUE)
# #   return(thisCountySum)
# # }
#
# #### Human footprint index - Venter ####
# # H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/HFP2009.tif
# hfiVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/HFP2009.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(hfiVenter)))
#
# hfiMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=hfiVenter)
# natureServeCountiesVenter$hfiMeans <- hfiMeans
# ## hfiSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = hfiVenter)
# ## natureServeCountiesVenter$hfiSums <- hfiSums
#
#
# #### Population density 1990 - Venter ####
# # H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Popdensity1990.tif
# popdens1990Venter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Popdensity1990.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(popdens1990Venter)))
#
# popdens1990Means <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=popdens1990Venter)
# natureServeCountiesVenter$popdens1990Means <- popdens1990Means
#
#
# #### Population density 2010 - Venter ####
# popdens2010Venter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Popdensity2010.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(popdens2010Venter)))
#
# popdens2010Means <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=popdens2010Venter)
# natureServeCountiesVenter$popdens2010Means <- popdens2010Means
#
#
# #### Built environments - Venter ####
# ## path: H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Built2009.tif
# builtVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Built1994.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(builtVenter)))
#
# builtMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=builtVenter)
# natureServeCountiesVenter$builtMeans <- builtMeans
# ## builtSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = builtVenter)
# ## natureServeCountiesVenter$builtSums <- builtSums
#
# #### Croplands - Venter ####
# ## path: H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/croplands2005.tif
# cropsVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/croplands2005.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(cropsVenter)))
#
# cropsMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=cropsVenter)
# natureServeCountiesVenter$cropsMeans <- cropsMeans
# ## cropsSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = cropsVenter)
# ## natureServeCountiesVenter$cropsSums <- cropsSums
#
# #### Pasture - Venter ####
# ## path: H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Pasture2009.tif
# pastureVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Pasture2009.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(pastureVenter)))
#
# pastureMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=pastureVenter)
# natureServeCountiesVenter$pastureMeans <- pastureMeans
# ## pastureSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = pastureVenter)
# ## natureServeCountiesVenter$pastureSums <- pastureSums
#
# #### Night lights - Venter ####
# ## path: H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Lights2009.tif
# lightVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Lights2009.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(lightVenter)))
#
# lightMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=lightVenter)
# natureServeCountiesVenter$nightLightMeans <- lightMeans
# ## lightSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = lightVenter)
# ## natureServeCountiesVenter$nightLightSums <- lightSums
#
# #### Roads - Venter ####
# ## path: H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Roads.tif
# roadsVenter <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/Maps/Roads.tif')
# if (!exists('natureServeCountiesVenter')) natureServeCountiesVenter <- spTransform(natureServeCounties, CRS(projection(roadsVenter)))
#
# roadMeans <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterMean, venterLayer=roadsVenter)
# natureServeCountiesVenter$roadMeans <- roadMeans
# ## lightSums <- sapply(natureServeCountiesVenter$FIPS_CODE_LONG, extractVenterSum, venterLayer = lightVenter)
# ## natureServeCountiesVenter$nightLightSums <- lightSums
#
#
# write.csv(natureServeCountiesVenter@data, "H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/nscWithVenterMeans.csv")
#
# rm(list=c('natureServeCountiesVenter', 'hifMeans', 'hfiVenter', 'builtMeans', 'builtVenter', 'popdens1990Venter', 'popdens1990Means', 'popdens2010Venter', 'popdens2010Means', 'cropsVenter', 'cropsMeans', 'pastureVenter', 'pastureMeans', 'lightVenter', 'lightMeans', 'roadsVenter', 'roadMeans','extractVenterMean'))


#### National land cover database ####
# nlcdRaster <- raster('H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img.img')
# natureServeCountiesNLCD <- spTransform(natureServeCounties, CRS(projection(nlcdRaster)))
#
# croppedNLCD <- crop(nlcdRaster, natureServeCountiesNLCD)
# rm(nlcdRaster)
#
# nlcdValue <- c(11,12,21,22,23,24,31,41,42,43,51,52,71,72,73,74,81,82,90,95)
# nlcdName <- c('open_water', 'ice_snow', 'developed_open','developed_low','developed_med', 'developed_high', 'barren','forest_deciduous','forest_evergreen', 'forest_mixed', 'scrub_dwarf','scrub_shrub','herb_grassland','herb_sedge','herb_lichens','herb_moss','cultivated_pasture','cultivated_crops', 'wetlands_woody','wetlands_emergent')
# nlcdLegend <- data.frame(nlcdValue, nlcdName)
#
# plot(croppedNLCD)
#
# nlcdTallies <- matrix(nrow  = 1, ncol = 21, data = c(0, nlcdValue))
#
# nrow(natureServeCountiesNLCD)
#
# for (i in 2001:3219) {
#   thisCounty <- natureServeCountiesNLCD[i, ]
#   if(thisCounty$STATE_NAME == 'Alaska') next
#   if (thisCounty$STATE_NAME == 'Hawaii') next
#   if (thisCounty$STATE_NAME == 'Puerto Rico') next
#   thisCountyCroppedNLCD <- crop(croppedNLCD, thisCounty)
#   plot(thisCountyCroppedNLCD)
#   thisCountyNLCD <- extract(thisCountyCroppedNLCD, thisCounty)
#   thisCountyNLCD <- unlist(thisCountyNLCD)
#
#   tallyNLCD <- function(x) {
#     cells <- length(thisCountyNLCD[ which(thisCountyNLCD == x) ])
#     #  cells <- rbind(x, cells)
#     return(cells)
#   }
#
#   thisCountyTallies <- sapply(nlcdLegend$nlcdValue, tallyNLCD)
#   thisCountyTallies <- unlist(thisCountyTallies)
#   thisCountyTallies <- c(as.numeric(thisCounty$FIPS_CODE_LONG), thisCountyTallies)
#   nlcdTallies <- rbind(nlcdTallies, thisCountyTallies)
#   write.csv(nlcdTallies, "H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies/nlcdTalliesSoFar2001.csv")
# }
#
# tallyFiles <- list.files(path = 'H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies', pattern = '.csv', all.files = TRUE, full.names = TRUE)
#
# for (i in 1:length(tallyFiles)){
#   if(!exists('nlcdAllTallies')) {
#     nlcdAllTallies <- read.csv(tallyFiles[1])
#     next
#   }
#   nlcdAllTallies <- rbind(nlcdAllTallies, read.csv(tallyFiles[i]))
#
# }
#
#
# colnames(nlcdAllTallies) <- nlcdAllTallies[1,]
# colnames(nlcdAllTallies) <- c('df', 'FIPS_CODE_LONG', colnames(nlcdAllTallies[3:length(nlcdAllTallies)]))
# nlcdAllTallies <- nlcdAllTallies[2:nrow(nlcdAllTallies), ]
# nlcdAllTallies <- unique(nlcdAllTallies)
#
# write.csv(nlcdAllTallies, 'H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies/nlcdAllTallies.csv')
# rm(nlcdAllTallies)
# rm(nlcdLegend)
# rm(nlcdTallies)
# rm(alaska)
# rm(croppedNLCD)
# rm(firststates)
# rm(hawaii)
# rm(n)
# rm(i)
# rm(nlcdName)
# rm(nlcdValue)
# rm(puertorico)
# rm(tallyFiles)
# rm(thisCountyNLCD)
# rm(thisCounty)
# rm(thisCountyTallies)
# rm(thisCountyCroppedNLCD)
# rm(tallyNLCD)
# rm(natureServeCountiesNLCD)
#

# load nlcd .csv
# nlcd <- read.csv('H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies/nlcdAllTallies.csv')
#
# countPixels <- function(i) {
#   nPixels <- sum(nlcd[i, 4:23])
#   return(nPixels)
# }
#
# nlcd$Total.Pixels <- sapply(1:nrow(nlcd), countPixels)
#
# getDevelopedFraction <- function(i) {
#   all.Dev.Sum <- sum(nlcd$X21[i], nlcd$X22[i], nlcd$X23[i])
#   all.Dev.Fraction <- all.Dev.Sum / nlcd$Total.Pixels[i]
#   return(all.Dev.Fraction)
# }
#
# nlcd$Developed.Fraction <- sapply(1:nrow(nlcd), getDevelopedFraction)
#
# getCultivatedFraction <- function(i) {
#   all.Cultivated.Sum <- sum(nlcd$X81[i], nlcd$X82[i])
#   all.Cultivated.Fraction <- all.Cultivated.Sum / nlcd$Total.Pixels[i]
#   return(all.Cultivated.Fraction)
# }
#
# nlcd$Cultivated.Fraction <- sapply(1:nrow(nlcd), getCultivatedFraction)
#
# getCropsFraction <- function(i) {
#   cropsFraction <- nlcd$X82[i]/nlcd$Total.Pixels[i]
#   return(cropsFraction)
# }
#
# nlcd$Crops.Fraction <- sapply(1:nrow(nlcd), getCropsFraction)
#
# getPastureFraction <- function(i) {
#   pastureFraction <- nlcd$X81[i]/nlcd$Total.Pixels[i]
#   return(pastureFraction)
# }
#
# nlcd$Pasture.Fraction <- sapply(1:nrow(nlcd), getPastureFraction)
#
# getGrasslandPastureFraction <- function(i) {
#   grass.Pasture.Sum <- sum(nlcd$X71[i], nlcd$X81[i])
#   grass.Pasture.Fraction <- grass.Pasture.Sum / nlcd$Total.Pixels[i]
# }
#
# nlcd$Grass.Pasture.Fraction <- sapply(1:nrow(nlcd), getGrasslandPastureFraction)
#
# write.csv(nlcd, 'H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies/nlcdAllTalliesWithFractions.csv')
# rm(countPixels)
# rm(getCropsFraction)
# rm(getCultivatedFraction)
# rm(getDevelopedFraction)
# rm(getGrasslandPastureFraction)
# rm(getPastureFraction)
# rm(nlcd)

# #### Mines ####
# activeMines <- shapefile('H:/Global Change Program/GIS/Land use/Mining/Active mines in the US/mineplant.shp')
# proj4string(activeMines)
# # natureServeCounties <- shapefile('H:/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/ORIGINAL/NS_mv_CTY_bdrys_G12ESAtots_201403')
# # proj4string(natureServeCounties)
#
# activeMines.2 <- spTransform(activeMines, CRS(proj4string(natureServeCounties)))
#
# mineCounties <- over(activeMines.2, natureServeCounties, returnList = FALSE)
# mineCounties <- cbind(mineCounties, activeMines.2@data)
#
# countMines <- function(FIPS_CODE_LONG) {
#   nmines <- nrow(mineCounties[ which(mineCounties$FIPS_CODE_LONG == FIPS_CODE_LONG), ])
#   return(nmines)
# }
#
# natureServeCounties_ActiveMines <- natureServeCounties
# natureServeCounties_ActiveMines$ACTIVE_MINES_COUNT <- sapply(natureServeCounties_ActiveMines$FIPS_CODE_LONG, countMines)
#
# ### These mines appear to be offshore
# mineCounties2 <- mineCounties[ which(is.na(mineCounties$FIPS_CODE_LONG)), ]
#
# activeMinesDf <- natureServeCounties_ActiveMines@data[ ,c('FIPS_CODE_LONG', "ACTIVE_MINES_COUNT")]
#
# write.csv(activeMinesDf, 'H:/Global Change Program/GIS/Land use/Mining/Active mines in the US/activeMines.csv')
# rm(activeMines)
# rm(activeMines.2)
# rm(countMines)
# rm(mineCounties)
# rm(mineCounties2)
# rm(activeMinesDf)
# rm(natureServeCounties_ActiveMines)




#### MRDS ####
### is another mines database ###

# natureServeCounties <- shapefile('H:/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/ORIGINAL/NS_mv_CTY_bdrys_G12ESAtots_201403')
# proj4string(natureServeCounties)
#
# mrds <- shapefile('H:/Global Change Program/GIS/Land use/Mining/mrds/mrds-trim.shp')
# proj4string(mrds)
#
# mrds.2 <- spTransform(mrds, CRS(proj4string(natureServeCounties)))
#
# mrds.2 <- mrds.2[ which(mrds.2$DEV_STAT == 'Producer' | mrds.2$DEV_STAT == 'Plant' | mrds.2$DEV_STAT == 'Past Producer'), ]
#
# mrdsCounties <- over(mrds.2, natureServeCounties, returnList = FALSE)
# mrdsCounties <- cbind(mrdsCounties, mrds.2@data)
#
# countMRDS <- function(FIPS_CODE_LONG) {
#   nmines <- nrow(mrdsCounties[ which(mrdsCounties$FIPS_CODE_LONG == FIPS_CODE_LONG), ])
#   return(as.numeric(nmines))
# }
#
# countedMines <- as.data.frame(natureServeCounties$FIPS_CODE_LONG)
# colnames(countedMines) <- c("FIPS_CODE_LONG")
# countedMines$MRDS_MINES_COUNT <- sapply(countedMines$FIPS_CODE_LONG, countMRDS)
#
# write.csv(countedMines, 'H:/Global Change Program/GIS/Land use/Mining/mrds/mrdsMinesProducers_Plants.csv')
#
# rm(mrds)
# rm(mrds.2)
# rm(mrdsCounties)
# rm(countMRDS)
# rm(countedMines)

#### Sriram et al layers ####
# Find original sources for these layers.

#### Srirarm cars ####
# path: H:/Global Change Program/GIS/Sriram/cars.tif
#
# sriramCars <- raster('H:/Global Change Program/GIS/Sriram/cars.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(sriramCars)))
# cars <- natureServeCountiesLayer$FIPS_CODE_LONG
# cars <- as.data.frame(cars)
# colnames(cars) <- 'FIPS_CODE_LONG'
# cars$meanCars <- NA
# for (i in 1:nrow(cars)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == cars$FIPS_CODE_LONG[i]), ]
#   cars$meanCars[i] <- extract(sriramCars, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# cars <- as.matrix(cars)
#
# write.csv(cars, 'H:/Global Change Program/GIS/Sriram/cars.csv')
#
# # most counties have 714.00 cars. Counties in Puerto Rico have -99 cars. May not be useful.
#
#
# #### Sriram soil degradation ####
# # path: H:/Global Change Program/GIS/Sriram/glasod.tif
# sriramSoilDeg <- raster('H:/Global Change Program/GIS/Sriram/glasod.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(sriramSoilDeg)))
# SoilDeg <- natureServeCountiesLayer$FIPS_CODE_LONG
# SoilDeg <- as.data.frame(SoilDeg)
# colnames(SoilDeg) <- 'FIPS_CODE_LONG'
# SoilDeg$meanSoilDeg <- NA
# for (i in 1:nrow(SoilDeg)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == SoilDeg$FIPS_CODE_LONG[i]), ]
#   SoilDeg$meanSoilDeg[i] <- extract(sriramSoilDeg, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# SoilDeg <- as.matrix(SoilDeg)
#
# write.csv(SoilDeg, 'H:/Global Change Program/GIS/Sriram/SoilDeg.csv')
#
#
# #### Sriram proximity to roads ####
# # path: H:/Global Change Program/GIS/Sriram/prox_roads.tif
# sriramProxRoads <- raster('H:/Global Change Program/GIS/Sriram/prox_roads.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(sriramProxRoads)))
# ProxRoads <- natureServeCountiesLayer$FIPS_CODE_LONG
# ProxRoads <- as.data.frame(ProxRoads)
# colnames(ProxRoads) <- 'FIPS_CODE_LONG'
# ProxRoads$meanProxRoads <- NA
# for (i in 1:nrow(ProxRoads)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == ProxRoads$FIPS_CODE_LONG[i]), ]
#   ProxRoads$meanProxRoads[i] <- extract(sriramProxRoads, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# ProxRoads <- as.matrix(ProxRoads)
#
# write.csv(ProxRoads, 'H:/Global Change Program/GIS/Sriram/ProxRoads.csv')
#
#
#
# #### Sriram pigs ####
# # path: H:/Global Change Program/GIS/Sriram/pigs.tif
# sriramPigs <- raster('H:/Global Change Program/GIS/Sriram/pigs.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(sriramPigs)))
# Pigs <- natureServeCountiesLayer$FIPS_CODE_LONG
# Pigs <- as.data.frame(Pigs)
# colnames(Pigs) <- 'FIPS_CODE_LONG'
# Pigs$meanPigs <- NA
# for (i in 1:nrow(Pigs)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == Pigs$FIPS_CODE_LONG[i]), ]
#   Pigs$meanPigs[i] <- extract(sriramPigs, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# Pigs <- as.matrix(Pigs)
#
# write.csv(Pigs, 'H:/Global Change Program/GIS/Sriram/Pigs.csv')
#
#
# #### Sriram poultry ####
# # path: H:/Global Change Program/GIS/Sriram/pigs.tif
# sriramPoultry <- raster('H:/Global Change Program/GIS/Sriram/poultry.tif')
# natureServeCountiesLayer <- spTransform(natureServeCounties, CRS(projection(sriramPoultry)))
# Poultry <- natureServeCountiesLayer$FIPS_CODE_LONG
# Poultry <- as.data.frame(Poultry)
# colnames(Poultry) <- 'FIPS_CODE_LONG'
# Poultry$meanPoultry <- NA
# for (i in 1:nrow(Poultry)) {
#   thisCounty <- natureServeCountiesLayer[ which(natureServeCountiesLayer$FIPS_CODE_LONG == Poultry$FIPS_CODE_LONG[i]), ]
#   Poultry$meanPoultry[i] <- extract(sriramPoultry, thisCounty, fun = mean, na.rm=TRUE)
# }
#
# Poultry <- as.matrix(Poultry)
#
# write.csv(Poultry, 'H:/Global Change Program/GIS/Sriram/Poultry.csv')

#### Zabel crop suitability ####
# path: H:/Global Change Program/GIS/Land use/Agricultural suitability/overall_cropsuit_i.shp
# path: H:/Global Change Program/GIS/Land use/Agricultural suitability/overall_cropsuit_i_1981-2010.tif

# cropSuit <- raster('H:/Global Change Program/GIS/Land use/Agricultural suitability/overall_cropsuit_i_1981-2010.tif')
# natureServeCountiesCropSuit <- spTransform(natureServeCounties, CRS(projection(cropSuit)))
# CropSuit <- natureServeCountiesCropSuit$FIPS_CODE_LONG
# CropSuit <- as.data.frame(CropSuit)
# colnames(CropSuit) <- 'FIPS_CODE_LONG'
# CropSuit$meanCropSuit <- NA
# for(i in 1:nrow(CropSuit)) {
#   thisCounty <- natureServeCountiesCropSuit[ which(natureServeCountiesCropSuit$FIPS_CODE_LONG == CropSuit$FIPS_CODE_LONG[i]), ]
#   CropSuit$meanCropSuit[i] <- extract(cropSuit, thisCounty, fun = mean, na.rm = TRUE)
# }
#
# CropSuit <- as.matrix(CropSuit)
# write.csv(CropSuit, 'H:/Global Change Program/GIS/Land use/Agricultural suitability/meanCropSuit.csv')
#
# rm(cropSuit)
# rm(CropSuit)
# rm(natureServeCountiesCropSuit)
# rm(thisCounty)
# rm(i)

#### Human appropriation of NPP ####
# path: H:/Global Change Program/GIS/Human Impact & Development/Human Appropriation of NPP - SEDAC/hapctnpp_geotiff.tif
#
# hanppRaster <- raster('H:/Global Change Program/GIS/Human Impact & Development/Human Appropriation of NPP - SEDAC/hapctnpp_geotiff.tif')
# plot(hanppRaster)
# natureServeCountiesHANPP <- spTransform(natureServeCounties, CRS(projection(hanppRaster)))
# hanppRaster <- crop(hanppRaster, natureServeCountiesHANPP)
# plot(hanppRaster)
#
# hanppDf <- natureServeCountiesHANPP$FIPS_CODE_LONG
# hanppDf <- as.data.frame(hanppDf)
# colnames(hanppDf) <- 'FIPS_CODE_LONG'
# hanppDf$meanHANPP <- NA
# for(i in 1:nrow(hanppDf)) {
#   thisCounty <- natureServeCountiesHANPP[ which(natureServeCountiesHANPP$FIPS_CODE_LONG == hanppDf$FIPS_CODE_LONG[i]), ]
#   hanppDf$meanHANPP[i] <- extract(hanppRaster, thisCounty, fun = mean, na.rm = TRUE)
# }
#
# hanppDf <- as.matrix(hanppDf)
# write.csv(hanppDf, 'H:/Global Change Program/GIS/Human Impact & Development/Human Appropriation of NPP - SEDAC/meanHANPP.csv')
#
# rm(hanppRaster)
# rm(natureServeCountiesHANPP)
# rm(hanppDf)
# rm(thisCounty)
# rm(i)

#### PADUS ####
# PADUS shapefile is too big to hold both the original shapefile and a subset in memory
#
# processPADUS_GAP <- function(gapLevel) {
#
#   padusShape <- shapefile('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/PADUS1_4Combined.shp')
#
#   padusShape <- padusShape[which(padusShape$GAP_Sts == as.character(gapLevel)), ]
#
#   plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
#   plainUSRaster <- plainUSRaster * 0
#
#   padusShape <- spTransform(padusShape, CRS(projection(plainUSRaster)))
#
#   plainUSRaster <- crop(plainUSRaster, padusShape)
#
#   padusRaster <- rasterize(padusShape, plainUSRaster)
#
#   padusRaster[!is.na(padusRaster)] <- 1
#
#   padusRasterMerged <- merge(padusRaster, plainUSRaster)
#
#   writeRaster(padusRasterMerged, paste0('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel', gapLevel), format ='GTiff', overwrite = TRUE, progress = 'text' )
# }
#
#
# processPADUS_GAP(1)
# pbPost("note", "PADUS complete", "GAP 1", recipients = c(1))
#
# processPADUS_GAP(2)
# pbPost("note", "PADUS complete", "GAP 2", recipients = c(1))
#
# processPADUS_GAP(3)
# pbPost("note", "PADUS complete", "GAP 3", recipients = c(1))
#
# processPADUS_GAP(4)
# pbPost("note", "PADUS complete", "GAP 4", recipients = c(1))
#
#
#
# processPADUS_ACC <- function(accessLevel) {
#
#   padusShape <- shapefile('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/PADUS1_4Combined.shp')
#
#   padusShape <- padusShape[which(padusShape$Access == as.character(accessLevel)), ]
#
#   plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
#   plainUSRaster <- plainUSRaster * 0
#
#   padusShape <- spTransform(padusShape, CRS(projection(plainUSRaster)))
#
#   plainUSRaster <- crop(plainUSRaster, padusShape)
#
#   padusRaster <- rasterize(padusShape, plainUSRaster)
#
#   padusRaster[!is.na(padusRaster)] <- 1
#
#   padusRasterMerged <- merge(padusRaster, plainUSRaster)
#
#   writeRaster(padusRasterMerged, paste0('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterAccessLevel', accessLevel), format ='GTiff', overwrite = TRUE, progress = 'text' )
# }
#
# processPADUS_ACC('UK')
# pbPost("note", "PADUS complete", "UK", recipients = c(1))
#
# processPADUS_ACC('RA')
# pbPost("note", "PADUS complete", "RA", recipients = c(1))
#
# processPADUS_ACC('OA')
# pbPost("note", "PADUS complete", "OA", recipients = c(1))
#
# processPADUS_ACC('XA')
# pbPost("note", "PADUS complete", "XA", recipients = c(1))
#
#
# processPADUS_IUCN <- function(iucn) {
#
#   padusShape <- shapefile('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/PADUS1_4Combined.shp')
#
#   padusShape <- padusShape[which(padusShape$IUCN_Cat == as.character(iucn)), ]
#
#   plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
#   plainUSRaster <- plainUSRaster * 0
#
#   padusShape <- spTransform(padusShape, CRS(projection(plainUSRaster)))
#
#   plainUSRaster <- crop(plainUSRaster, padusShape)
#
#   padusRaster <- rasterize(padusShape, plainUSRaster)
#
#   padusRaster[!is.na(padusRaster)] <- 1
#
#   padusRasterMerged <- merge(padusRaster, plainUSRaster)
#
#   writeRaster(padusRasterMerged, paste0('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevel', iucn), format ='GTiff', overwrite = TRUE, progress = 'text' )
# }
#
# processPADUS_IUCN('Other Conservation Area')
# pbPost("note", "PADUS complete", "OCA", recipients = c(1))
#
# processPADUS_IUCN('IV')
# pbPost("note", "PADUS complete", "IV", recipients = c(1))
#
# processPADUS_IUCN('Ia')
# pbPost("note", "PADUS complete", "Ia", recipients = c(1))
#
# processPADUS_IUCN('V')
# pbPost("note", "PADUS complete", "V" recipients = c(1))
#
# processPADUS_IUCN('VI')
# pbPost("note", "PADUS complete", "VI", recipients = c(1))
#
# processPADUS_IUCN('Unassigned')
# pbPost("note", "PADUS complete", "Unassigned", recipients = c(1))
#
# processPADUS_IUCN('III')
# pbPost("note", "PADUS complete", "III", recipients = c(1))
#
# processPADUS_IUCN('II')
# pbPost("note", "PADUS complete", "II", recipients = c(1))
#
# processPADUS_IUCN('Ib')
# pbPost("note", "PADUS complete", "Ib", recipients = c(1))
#
# processPADUS_IUCN('N/R')
# pbPost("note", "PADUS complete", "NR", recipients = c(1))
#
#
#
# # # end PADUS shapefile processing
# # from now on use rasters
# padusGAPLevel1 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel1.tif')
#
# natureServeCountiesPADUS <- spTransform(natureServeCounties, CRS(projection(padusGAPLevel1)))
#
# padusProportions <- as.data.frame(natureServeCountiesPADUS$FIPS_CODE_LONG)
# colnames(padusProportions) <- c('FIPS_CODE_LONG')
#
#
# proportionPA <- function(FIPS_CODE_LONG, rasterObject) {
#   thisCounty <- natureServeCountiesPADUS[ which(natureServeCountiesPADUS$FIPS_CODE_LONG == FIPS_CODE_LONG), ]
#   rasterObject <- crop(rasterObject, extent(thisCounty))
#   thisCountyProportion <- extract(rasterObject, thisCounty, method='simple', fun = mean, na.rm = TRUE)
#   if (is.null(thisCountyProportion[[1]])) thisCountyProportion <- NA
#   if (!is.null(thisCountyProportion[[1]])) thisCountyProportion <- as.numeric(thisCountyProportion)
#   #}
#   return(thisCountyProportion)
# }
#
#
# gap1Proportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusGAPLevel1)
# gap1Proportion <- unlist(gap1Proportion)
# padusProportions$gap1Proportion <- gap1Proportion
#
# padusGAPLevel2 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel2.tif')
# gap2Proportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusGAPLevel2)
# gap2Proportion <- unlist(gap2Proportion)
# padusProportions$gap2Proportion <- gap2Proportion
#
# padusGAPLevel3 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel3.tif')
# gap3Proportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusGAPLevel3)
# gap3Proportion <- unlist(gap3Proportion)
# padusProportions$gap3Proportion <- gap3Proportion
#
# padusGAPLevel4 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel4.tif')
# gap4Proportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusGAPLevel4)
# gap4Proportion <- unlist(gap4Proportion)
# padusProportions$gap4Proportion <- gap4Proportion
#
#
# rm(padusGAPLevel1)
# rm(padusGAPLevel2)
# rm(padusGAPLevel3)
# rm(padusGAPLevel4)
# rm(gap2Proportion)
# rm(gap3Proportion)
# rm(gap4Proportion)
# rm(gap1Proportion)
# pbPost("note", "PADUS complete", "GAP extractions", recipients = c(1))
#
# padusAccessUK <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterAccessLevelUK.tif')
# accessUKProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusAccessUK)
# accessUKProportion <- unlist(accessUKProportion)
# padusProportions$accessUKProportion <- accessUKProportion
# rm(accessUKProportion)
#
# padusAccessRA <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterAccessLevelRA.tif')
# accessRAProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusAccessRA)
# accessRAProportion <- unlist(accessRAProportion)
# padusProportions$accessRAProportion <- accessRAProportion
# rm(accessRAProportion)
#
#
# padusAccessOA <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterAccessLevelOA.tif')
# accessOAProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusAccessOA)
# accessOAProportion <- unlist(accessOAProportion)
# padusProportions$accessOAProportion <- accessOAProportion
# rm(accessOAProportion)
#
#
# padusAccessXA <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterAccessLevelXA.tif')
# accessXAProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusAccessXA)
# accessXAProportion <- unlist(accessXAProportion)
# padusProportions$accessXAProportion <- accessXAProportion
# rm(accessXAProportion)
#
# pbPost("note", "PADUS complete", "Access extractions", recipients = c(1))
#
#
# padusIUCN_Ia <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelIa.tif')
# IUCN_IaProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_Ia)
# IUCN_IaProportion <- unlist(IUCN_IaProportion)
# padusProportions$IUCN_IaProportion <- IUCN_IaProportion
# rm(IUCN_IaProportion)
#
#
# padusIUCN_Ib <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelIb.tif')
# IUCN_IbProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_Ib)
# IUCN_IbProportion <- unlist(IUCN_IbProportion)
# padusProportions$IUCN_IbProportion <- IUCN_IbProportion
# rm(IUCN_IbProportion)
#
#
# padusIUCN_II <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelII.tif')
# IUCN_IIProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_II)
# IUCN_IIProportion <- unlist(IUCN_IIProportion)
# padusProportions$IUCN_IIProportion <- IUCN_IIProportion
# rm(IUCN_IIProportion)
#
#
# padusIUCN_III <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelIII.tif')
# IUCN_IIIProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_III)
# IUCN_IIIProportion <- unlist(IUCN_IIIProportion)
# padusProportions$IUCN_IIIProportion <- IUCN_IIIProportion
# rm(IUCN_IIIProportion)
#
#
# padusIUCN_IV <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelIV.tif')
# IUCN_IVProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_IV)
# IUCN_IVProportion <- unlist(IUCN_IVProportion)
# padusProportions$IUCN_IVProportion <- IUCN_IVProportion
# rm(IUCN_IVProportion)
#
#
# padusIUCN_V <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelV.tif')
# IUCN_VProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_V)
# IUCN_VProportion <- unlist(IUCN_VProportion)
# padusProportions$IUCN_VProportion <- IUCN_VProportion
# rm(IUCN_VProportion)
#
#
# padusIUCN_VI <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelVI.tif')
# IUCN_VIProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_VI)
# IUCN_VIProportion <- unlist(IUCN_VIProportion)
# padusProportions$IUCN_VIProportion <- IUCN_VIProportion
# rm(IUCN_VIProportion)
#
#
# padusIUCN_Other <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelOther Conservation Area.tif')
# IUCN_OtherProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_Other)
# IUCN_OtherProportion <- unlist(IUCN_OtherProportion)
# padusProportions$IUCN_OtherProportion <- IUCN_OtherProportion
# rm(IUCN_OtherProportion)
#
#
# padusIUCN_Unassgn <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelUnassigned.tif')
# IUCN_UnassgnProportion <- lapply(natureServeCountiesPADUS$FIPS_CODE_LONG, proportionPA, rasterObject = padusIUCN_Unassgn)
# IUCN_UnassgnProportion <- unlist(IUCN_UnassgnProportion)
# padusProportions$IUCN_UnassgnProportion <- IUCN_UnassgnProportion
# rm(IUCN_UnassgnProportion)
#
# pbPost("note", "PADUS complete", "IUCN extractions", recipients = c(1))
#
# write.csv(padusProportions, "H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/PADUS_Proportions.csv")
#
# # Get target columns
#
# # GAP 1 + 2  + 3
#
# gap1 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel1.tif')
# gap2 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel2.tif')
# gap3 <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterGAPLevel3.tif')
# iucn1a <- raster('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/Rasters/padusRasterIUCNLevelIa.tif')
#
#
# natureServeCountiesPADUS <- spTransform(natureServeCounties, CRS(projection(gap1)))
#
# getPADUSMeans <- function(FIPS_CODE_LONG) {
#
#   thisCounty <- natureServeCountiesPADUS[ which(natureServeCountiesPADUS$FIPS_CODE_LONG == FIPS_CODE_LONG), ]
#
#   this.g1 <- crop(gap1, extent(thisCounty))
#   this.g1 <- extract(this.g1, thisCounty, method = 'simple')
#   this.g1 <- as.vector(unlist(this.g1))
#   this.g2 <- crop(gap2, extent(thisCounty))
#   this.g2 <- extract(this.g2, thisCounty, method = 'simple')
#   this.g2 <- as.vector(unlist(this.g2))
#   this.g3 <- crop(gap3, extent(thisCounty))
#   this.g3 <- extract(this.g3, thisCounty, method = 'simple')
#   this.g3 <- as.vector(unlist(this.g3))
#   this.1a <- crop(iucn1a, extent(thisCounty))
#   this.1a <- extract(this.1a, thisCounty, method = 'simple')
#   this.1a <- as.vector(unlist(this.1a))
#
#   this.g1 <- as.logical(this.g1)
#   this.g2 <- as.logical(this.g2)
#   this.g3 <- as.logical(this.g3)
#   this.1a <- as.logical(this.1a)
#
#   this.12 <- this.g1 | this.g2
#   this.123 <- this.g1 | this.g2 | this.g3
#   this.3X12 <- this.g3 & !this.12
#   this.12X1a <- this.12 & !this.1a
#   this.123X1a <- this.123 & !this.1a
#
#   this <- cbind(this.12, this.123, this.3X12, this.12X1a, this.123X1a, this.1a)
#   this <- as.matrix(this)
#   this <- as.data.frame(this)
#
#   this.12 <- mean(this$this.12, na.rm = TRUE)
#   this.123 <- mean(this$this.123, na.rm = TRUE)
#   this.3X12 <- mean(this$this.3X12, na.rm = TRUE)
#   this.12X1a <- mean(this$this.12X1a, na.rm = TRUE)
#   this.123X1a <- mean(this$this.123X1a, na.rm = TRUE)
#   this.1a <- mean(this$this.1a, na.rm = TRUE)
#
#   out <- list (iucn1a = this.1a,
#                gap12 = this.12,
#                gap123 = this.123X1a,
#                gap3X12 = this.3X12,
#                gap12X1a = this.12X1a,
#                gap123X1a = this.123X1a)
#   return(out)
# }
#
# outs <- list(length=nrow(natureServeCountiesPADUS@data))
#
# for(i in 1:nrow(natureServeCountiesPADUS@data)) {
#   outs[[i]] <- getPADUSMeans(natureServeCountiesPADUS$FIPS_CODE_LONG[i])
# }
#
# padus <- matrix(data = NA, nrow = (nrow(natureServeCountiesPADUS@data)), ncol = 7)
# padus <- as.data.frame(padus)
# colnames(padus) <- c('FIPS_CODE_LONG', 'iucn1a', 'gap12', 'gap123', 'gap3X12', 'gap12X1a', 'gap123X1a')
#
# for(i in 1:nrow(natureServeCountiesPADUS@data)) {
#   padus[i, 1] <- natureServeCountiesPADUS$FIPS_CODE_LONG[i]
#   padus[i, 2:7] <- as.vector(unlist(outs[i]))
# }
#
# write.csv(padus, 'H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/padusCounties.csv')


#### Invasive alien species - EDDMapS ####
#
# Don't use for now - permissions unclear
# iasFiles <- list.files(path = 'H:/Global Change Program/GIS/Invasive species/EDDMapS April 12 2017', full.names = TRUE)
#
# readIAS <- function(iasFilePath) {
#   tester <- read.table(iasFilePath, nrows = 1)
#   if (tester[1]!="Parameters") return()
#   iasHeader <- read.csv(iasFilePath, nrows = 3)
#   iasData <- read.csv(iasFilePath, skip = 5)
#   colnames(iasData)[2] <- paste0('Number.of.Invasive.', iasHeader$Parameters[2])
#   iasData$State <- iasHeader$Parameters[1]
#   return(iasData)
# }
#
# allIAS <- lapply(iasFiles, readIAS)
#
#
# masterDf <- natureServeCounties@data[ , c('FIPS_CODE_LONG', 'NAME', 'STATE_NAME')]
# masterDf$Number.of.Invasive.Diseases <- NA
# masterDf$Number.of.Invasive.Plants <- NA
# masterDf$Number.of.Invasive.Wildlife <- NA
# masterDf$Number.of.Invasive.Insects <- NA
#
#
#
# combinedIAS <- as.data.frame(allIAS[1])
# combinedIAS$Number.of.Invasive.Diseases <- NA
# combinedIAS$Number.of.Invasive.Plants <- NA
# combinedIAS$Number.of.Invasive.Wildlife <- NA
# combinedIAS <- combinedIAS[0, ]
#
# for (i in 1:length(allIAS)) {
#   thisDf <- as.data.frame(allIAS[i])
#   if(nrow(thisDf) < 1) next
#   if(colnames(thisDf)[2] == 'Number.of.Invasive.Insects') {
#     thisDf$Number.of.Invasive.Diseases <- NA
#     thisDf$Number.of.Invasive.Plants <- NA
#     thisDf$Number.of.Invasive.Wildlife <- NA
#   }
#
#   if(colnames(thisDf)[2] == 'Number.of.Invasive.Diseases') {
#     thisDf$Number.of.Invasive.Insects <- NA
#     thisDf$Number.of.Invasive.Plants <- NA
#     thisDf$Number.of.Invasive.Wildlife <- NA
#   }
#
#   if(colnames(thisDf)[2] == 'Number.of.Invasive.Plants') {
#     thisDf$Number.of.Invasive.Diseases <- NA
#     thisDf$Number.of.Invasive.Insects <- NA
#     thisDf$Number.of.Invasive.Wildlife <- NA
#   }
#
#   if(colnames(thisDf)[2] == 'Number.of.Invasive.Wildlife') {
#     thisDf$Number.of.Invasive.Diseases <- NA
#     thisDf$Number.of.Invasive.Plants <- NA
#     thisDf$Number.of.Invasive.Insects <- NA
#   }
#
#   if(exists('combinedIAS')) combinedIAS <- rbind(combinedIAS, thisDf)
#   if(!exists('combinedIAS')) combinedIAS <- thisDf
# }
#
# ### for every line in combinedIAS...
#
#
# combinedIAS$FIPS_CODE_LONG <- NA
# combinedIAS$rowIndex <- c(1:nrow(combinedIAS))
# problemCounties <- vector()
#
# for (i in 1:nrow(combinedIAS)){
#   # figure out the corresponding row and figure out the corresponding column...
#
#   justCty <- unlist(strsplit(combinedIAS$County.Name[i], split = " County"))[1]
#   justCty <- unlist(strsplit(justCty, split = " city"))[1]
#   justCty <- unlist(strsplit(justCty, split = " Cens"))[1]
#   justCty <- unlist(strsplit(justCty, split = " Pari"))[1]
#   justCty <- unlist(strsplit(justCty, split = " Boro"))[1]
#   justCty <- unlist(strsplit(justCty, split = " City"))[1]
#   if (justCty == 'Lagrange') justCty <- 'LaGrange'
#   if (justCty == 'Hoonah-Angoon') justCty <- 'Skagway-Hoonah-Angoon'
#   if (justCty == 'Skagway Municipality') justCty <- 'Skagway-Hoonah-Angoon'
#   if (justCty == 'Prince of Wales-Outer Ket') justCty <- 'Prince of Wales-Outer Ketchikan'
#   if (justCty == 'Wrangell') justCty <- 'Wrangell-Petersburg'
#   if (justCty == 'De Kalb') justCty <- 'DeKalb'
#   if (justCty == 'La Porte') justCty <- 'LaPorte'
#   if (justCty == 'DeBaca') justCty <- 'De Baca'
#   if (justCty == 'Mc Kean') justCty <- 'McKean'
#   if (justCty == 'James') justCty <- 'James City'
#
#
#
#   ctyOptions <- masterDf[which(masterDf$NAME == justCty), ]
#   if (nrow(ctyOptions) > 1){
#     ctyOptions <- ctyOptions[which(ctyOptions$STATE_NAME == combinedIAS$State[i]), ]
#   }
#   if (nrow(ctyOptions) < 1){
#     problemCounties <- c(problemCounties, i)
#     next
#   }
#
#   combinedIAS$FIPS_CODE_LONG[i] <- ctyOptions$FIPS_CODE_LONG
# }
#
# # Prince of Wales - Hyder Census Area is not in NSC.
# # Escambia Cty, Illinois is not in NSC (and does not exist?)
# # Washington Cty, Washington is not in NSC (and does not exist?)
#
#
#
# invasivePlants <- combinedIAS[ , c('rowIndex', 'County.Name', 'State', 'Number.of.Invasive.Plants')]
# invasivePlants <- na.omit(invasivePlants)
# invasivePlants <- join(invasivePlants, combinedIAS[ ,c('rowIndex', 'FIPS_CODE_LONG')], by = 'rowIndex', type = 'left')
#
#
#
# invasiveDiseases <- combinedIAS[ , c('rowIndex', 'County.Name', 'State', 'Number.of.Invasive.Diseases')]
# invasiveDiseases <- na.omit(invasiveDiseases)
# invasiveDiseases <- join(invasiveDiseases, combinedIAS[ ,c('rowIndex', 'FIPS_CODE_LONG')], by = 'rowIndex', type = 'left')
#
# invasiveInsects <- combinedIAS[ , c('rowIndex', 'County.Name', 'State', 'Number.of.Invasive.Insects')]
# invasiveInsects <- na.omit(invasiveInsects)
# invasiveInsects <- join(invasiveInsects, combinedIAS[ ,c('rowIndex', 'FIPS_CODE_LONG')], by = 'rowIndex', type = 'left')
#
# invasiveWildlife <- combinedIAS[ , c('rowIndex', 'County.Name', 'State', 'Number.of.Invasive.Wildlife')]
# invasiveWildlife <- na.omit(invasiveWildlife)
# invasiveWildlife <- join(invasiveWildlife, combinedIAS[ ,c('rowIndex', 'FIPS_CODE_LONG')], by = 'rowIndex', type = 'left')
#
# natureServeCountiesInvasives <- natureServeCounties
# natureServeCountiesInvasives@data <- join(natureServeCountiesInvasives@data, invasivePlants[,c('FIPS_CODE_LONG', 'Number.of.Invasive.Plants')], by = 'FIPS_CODE_LONG', type = 'left')
# natureServeCountiesInvasives@data <- join(natureServeCountiesInvasives@data, invasiveDiseases[,c('FIPS_CODE_LONG', 'Number.of.Invasive.Diseases')], by = 'FIPS_CODE_LONG', type = 'left')
# natureServeCountiesInvasives@data <- join(natureServeCountiesInvasives@data, invasiveInsects[,c('FIPS_CODE_LONG', 'Number.of.Invasive.Insects')], by = 'FIPS_CODE_LONG', type = 'left')
# natureServeCountiesInvasives@data <- join(natureServeCountiesInvasives@data, invasiveWildlife[,c('FIPS_CODE_LONG', 'Number.of.Invasive.Wildlife')], by = 'FIPS_CODE_LONG', type = 'left')
#
# write.csv(natureServeCountiesInvasives@data, 'H:/Global Change Program/GIS/Invasive species/EDDMapS April 12 2017/nsvWithInvasives.csv')


#### NLCD Percent Developed Imperviousness ####
# nlcdImpervious <- raster('H:/Global Change Program/GIS/Land use/NCLD Developed Imperviousness 2011/impe48i0100a.tif')
# natureServeCountiesImpervious <- spTransform(natureServeCounties, CRS(projection(nlcdImpervious)))
# natureServeCountiesImpervious$Percent.Impervious <- NA
#
#
# for (i in 1:nrow(natureServeCountiesImpervious)) {
#   percent <- extract(nlcdImpervious, natureServeCountiesImpervious[i, ], fun = mean, na.omit =  TRUE)
#   if (is.null(percent[[1]])) percent <- NA
#   natureServeCountiesImpervious$Percent.Impervious[i] <- percent
#   write.csv(natureServeCountiesImpervious@data, 'H:/Global Change Program/GIS/Land use/NCLD Developed Imperviousness 2011/percentImpervious.csv')
#
# }
#
# write.csv(natureServeCountiesImpervious@data, 'H:/Global Change Program/GIS/Land use/NCLD Developed Imperviousness 2011/percentImpervious.csv')
#

# #### SPARROW Water quality (N and P pollution) ####
#
# # New England and Lower Atlantic use one set of shapefiles...
# # The catches span 2 shapefiles. They are too large to load all at once, so I am rasterizing each one separately.
#
# nitrogenNELA <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 NE and Mid Atlantic.csv', skip = 121)
# nitrogenNELA <- nitrogenNELA[, 1:(ncol(nitrogenNELA) - 1)]
# colnames(nitrogenNELA) <- c('COMID', colnames(nitrogenNELA)[2:length(colnames(nitrogenNELA))])
#
# nitrogenNELA <- nitrogenNELA[ , c(1, 2, 3, 5, 12, 19, 20, 21)]
# colnames(nitrogenNELA) <- c('COMID', 'Total_Contributing_Area', 'Total_Upstream_Area', "Total_Yield", 'Incremental_Load', 'Total_Load', 'HUC8', 'Reach_Name')
#
# catch1 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/Catchments/April25/NHDPlus01/Drainage/catchment.shp')
# catch1@data <- join(catch1@data, nitrogenNELA, by = 'COMID', type = 'left')
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch1 <- spTransform(catch1, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch1)
#
# totalYield1.r <- rasterize(catch1, plainUSRaster, "Total_Yield")
# # incYield1.r <- rasterize(catch1, plainUSRaster, "Incremental_Yield")
#
# writeRaster(totalYield1.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen rasters/totalYield1', format ='GTiff', overwrite = TRUE, progress = 'text' )
# # writeRaster(incYield1.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen rasters/incLoad1', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(totalYield1.r)
# # rm(incLoad1.r)
# rm(plainUSRaster)
# rm(catch1)
# rm(keep)
#
# catch2 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/Catchments/April25/NHDPlus02/Drainage/catchment.shp')
# catch2@data <- join(catch2@data, nitrogenNELA, by = 'COMID', type = 'left')
# keep <- which(!is.na(catch2$Total_Yield))
# catch2 <- catch2 [keep, ]
#
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch2 <- spTransform(catch2, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch2)
#
# totalYield2.r <- rasterize(catch2, plainUSRaster, "Total_Yield")
# # incYield2.r <- rasterize(catch2, plainUSRaster, "Incremental_Yield")
#
# writeRaster(totalYield2.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield2', format ='GTiff', overwrite = TRUE, progress = 'text' )
# # writeRaster(incYield2.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/incYield2', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(totalYield2.r)
# # rm(incYield2.r)
# rm(plainUSRaster)
# rm(catch2)
# rm(keep)
# rm(nitrogenNELA)
#
# ### All other river basins use the same shapefile and can be handled simultaneously
#
# # lower Mississippi basin
# nitrogenLowMiss <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 Lower Mississippi etc.csv', skip = 121)
# colnames(nitrogenLowMiss) <- c('MRB_ID', colnames(nitrogenLowMiss)[2:length(colnames(nitrogenLowMiss))])
# nitrogenLowMiss <- nitrogenLowMiss[,1:(ncol(nitrogenLowMiss) -1)]
#
# nitrogenLowMiss <- nitrogenLowMiss[ , c(1, 2, 3, 5, 13, 21, 22, 23)]
#
# # Missouri basin
# nitrogenMissouri <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 Missouri.csv', skip = 121)
# colnames(nitrogenMissouri) <- c('MRB_ID', colnames(nitrogenMissouri)[2:length(colnames(nitrogenMissouri))])
# nitrogenMissouri <- nitrogenMissouri[,1:(ncol(nitrogenMissouri) - 1)]
#
# nitrogenMissouri <- nitrogenMissouri[ , colnames(nitrogenLowMiss)]
#
# # Pacific Northwest
# nitrogenPNW <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 Pacific Northwest.csv', skip = 121)
# colnames(nitrogenPNW) <- c('MRB_ID', colnames(nitrogenPNW)[2:length(colnames(nitrogenPNW))])
# nitrogenPNW <- nitrogenPNW[ , 1:(ncol(nitrogenPNW) - 1)]
#
# nitrogenPNW <- nitrogenPNW[ , colnames(nitrogenLowMiss)]
#
#
# # South Atlantic Gulf
# nitrogenSAG <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 South Atlantic Gulf and Tennessee.csv', skip = 121)
# colnames(nitrogenSAG) <- c('MRB_ID', colnames(nitrogenSAG)[2:length(colnames(nitrogenSAG))])
# nitrogenSAG <- nitrogenSAG[ , 1:(ncol(nitrogenSAG) - 1)]
#
# nitrogenSAG <- nitrogenSAG[ , colnames(nitrogenLowMiss)]
#
# # Upper Midwest
# nitrogenUMW <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Nitrogen 2002 Upper Midwest.csv', skip = 121)
# colnames(nitrogenUMW) <- c('MRB_ID', colnames(nitrogenUMW)[2:length(colnames(nitrogenUMW))])
# nitrogenUMW <- nitrogenUMW[ , 1:(ncol(nitrogenUMW) - 1)]
#
# nitrogenUMW <- nitrogenUMW[ , colnames(nitrogenLowMiss)]
#
#
#
# nitrogenAll <- rbind(nitrogenMissouri, nitrogenPNW, nitrogenSAG, nitrogenUMW, nitrogenLowMiss)
#
# rm(nitrogenMissouri)
# rm(nitrogenPNW)
# rm(nitrogenSAG)
# rm(nitrogenUMW)
# rm(nitrogenLowMiss)
#
# colnames(nitrogenAll) <- c("MRB_ID", 'Total_Contributing_Area', 'Total_Upstream_Area', 'Total_Yield', 'Incremental_Load', 'Total_Load', 'HUC8', 'Reach_Name')
#
# # No NA's in nitrogenAll
# anyNA(nitrogenAll)
#
# catch3 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/2002_sparrow_MRB23457.shp')
#
# # Catches with no corresponding shapefile also have no HUC8 code. There are 12 of them.
# catch3.test <- catch3
# catch3.test@data$Shape <- 1
# catchesNoShape <- join(nitrogenAll, catch3.test@data, by = "MRB_ID", type = 'left')
# nrow(catchesNoShape[ which(is.na(catchesNoShape$Shape)), ])
# rm(catchesNoShape)
# rm(catch3.test)
#
# # One catch in nitrogenAll has two records
# # Remove the second, smaller, record
# notUnique <- vector()
# MRBS <- nitrogenAll$MRB_ID
#
# for (i in 1:length(MRBS)) {
#   thisMRB <- MRBS[i]
#   matches <- which(MRBS == thisMRB)
#   if(length(matches) > 1) notUnique <- c(notUnique, i)
#
#
# }
#
# nitrogenAll[32654, 'Total_Yield']
# nitrogenAll[37851, ]
#
# nitrogenAll <- nitrogenAll[c(1:37850, 37852:nrow(nitrogenAll)), ]
#
#
# # Join catches
# catch3@data <- join(catch3@data, nitrogenAll, by = "MRB_ID", type = 'left')
# # There are no catches without modeled values.
# anyNA(catch3$Total_Yield)
#
#
#
# View(catch3@data)
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch3 <- spTransform(catch3, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch3)
#
# totalYield3.r <- rasterize(catch3, plainUSRaster, "Total_Yield")
# # incLoad3.r <- rasterize(catch3, plainUSRaster, "Incremental_Yield")
#
# writeRaster(totalYield3.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield3', format ='GTiff', overwrite = TRUE, progress = 'text' )
# # writeRaster(incLoad3.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/incLoad3', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(nitrogenAll)
# rm(catch3)
# rm(plainUSRaster)
# rm(totalYield3.r)
#
#
# nitrogen1 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield1.tif')
# nitrogen2 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield2.tif')
# nitrogen3 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield3.tif')
# ?raster
#
# plot(intersect(nitrogen1, nitrogen2))
#
# nitrogenTogether <- mosaic(nitrogen1, nitrogen2, fun = max)
# nitrogenTogether <- mosaic(nitrogenTogether, nitrogen3, fun = max)
# writeRaster(nitrogenTogether,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield_US', format ='GTiff', overwrite = TRUE, progress = 'text' )
# rm(nitrogen1)
# rm(nitrogen2)
# rm(nitrogen3)
# rm(nitrogenTogether)
#
# ### Phosphorus ###
#
# # New England and Lower Atlantic use one set of shapefiles...
# # The catches span 2 shapefiles. They are too large to load all at once, so I am rasterizing each one separately.
#
# phosphorusNELA <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus NE and Mid Atlantic.csv', skip = 121)
# phosphorusNELA <- phosphorusNELA[, 1:(ncol(phosphorusNELA) - 1)]
# colnames(phosphorusNELA) <- c('COMID', colnames(phosphorusNELA)[2:length(colnames(phosphorusNELA))])
#
# phosphorusNELA <- phosphorusNELA[ , c(1, 2, 3, 5, 12, 19, 20, 21)]
# colnames(phosphorusNELA)
# colnames(phosphorusNELA) <- c('COMID', 'Total_Contributing_Area', 'Total_Upstream_Area', "Total_Yield", 'Incremental_Load', 'Total_Load', 'HUC8', 'Reach_Name')
#
# catch1 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/Catchments/April25/NHDPlus01/Drainage/catchment.shp')
# catch1@data <- join(catch1@data, phosphorusNELA, by = 'COMID', type = 'left')
#
#
# anyNA(catch1@data$Total_Yield)
#
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch1 <- spTransform(catch1, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch1)
#
# totalYield1.r <- rasterize(catch1, plainUSRaster, "Total_Yield")
# # incYield1.r <- rasterize(catch1, plainUSRaster, "Incremental_Yield")
#
# writeRaster(totalYield1.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus rasters/totalYield1', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(totalYield1.r)
# # rm(incLoad1.r)
# rm(plainUSRaster)
# rm(catch1)
# rm(keep)
#
# catch2 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/Catchments/April25/NHDPlus02/Drainage/catchment.shp')
# catch2@data <- join(catch2@data, phosphorusNELA, by = 'COMID', type = 'left')
# keep <- which(!is.na(catch2$Total_Yield))
# catch2 <- catch2 [keep, ]
#
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch2 <- spTransform(catch2, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch2)
#
# totalYield2.r <- rasterize(catch2, plainUSRaster, "Total_Yield")
# # incYield2.r <- rasterize(catch2, plainUSRaster, "Incremental_Yield")
#
# plot(totalYield2.r)
#
# writeRaster(totalYield2.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield2', format ='GTiff', overwrite = TRUE, progress = 'text' )
# # writeRaster(incYield2.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/incYield2', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(totalYield2.r)
# # rm(incYield2.r)
# rm(plainUSRaster)
# rm(catch2)
# rm(keep)
# rm(phosphorusNELA)
#
# ### All other river basins use the same shapefile and can be handled simultaneously
#
# # lower Mississippi basin
# phosphorusLowMiss <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus Lower Mississippi.csv', skip = 121)
# colnames(phosphorusLowMiss) <- c('MRB_ID', colnames(phosphorusLowMiss)[2:length(colnames(phosphorusLowMiss))])
# phosphorusLowMiss <- phosphorusLowMiss[,1:(ncol(phosphorusLowMiss) -1)]
#
# phosphorusLowMiss <- phosphorusLowMiss[ , c(1, 2, 3, 5, 13, 21, 22, 23)]
#
# colnames(phosphorusLowMiss)
#
# # Missouri basin
# phosphorusMissouri <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus Missouri.csv', skip = 121)
# colnames(phosphorusMissouri) <- c('MRB_ID', colnames(phosphorusMissouri)[2:length(colnames(phosphorusMissouri))])
# phosphorusMissouri <- phosphorusMissouri[,1:(ncol(phosphorusMissouri) - 1)]
#
# phosphorusMissouri <- phosphorusMissouri[ , colnames(phosphorusLowMiss)]
#
# # Pacific Northwest
# phosphorusPNW <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus Pacific Northwest.csv', skip = 121)
# colnames(phosphorusPNW) <- c('MRB_ID', colnames(phosphorusPNW)[2:length(colnames(phosphorusPNW))])
# phosphorusPNW <- phosphorusPNW[ , 1:(ncol(phosphorusPNW) - 1)]
#
# phosphorusPNW <- phosphorusPNW[ , colnames(phosphorusLowMiss)]
#
#
# # South Atlantic Gulf
# phosphorusSAG <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus Southeast.csv', skip = 121)
# colnames(phosphorusSAG) <- c('MRB_ID', colnames(phosphorusSAG)[2:length(colnames(phosphorusSAG))])
# phosphorusSAG <- phosphorusSAG[ , 1:(ncol(phosphorusSAG) - 1)]
#
# phosphorusSAG <- phosphorusSAG[ , colnames(phosphorusLowMiss)]
#
# # Upper Midwest
# phosphorusUMW <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Total Phosphorus Upper Midwest.csv', skip = 121)
# colnames(phosphorusUMW) <- c('MRB_ID', colnames(phosphorusUMW)[2:length(colnames(phosphorusUMW))])
# phosphorusUMW <- phosphorusUMW[ , 1:(ncol(phosphorusUMW) - 1)]
#
# phosphorusUMW <- phosphorusUMW[ , colnames(phosphorusLowMiss)]
#
#
#
# phosphorusAll <- rbind(phosphorusMissouri, phosphorusPNW, phosphorusSAG, phosphorusUMW, phosphorusLowMiss)
#
# rm(phosphorusMissouri)
# rm(phosphorusPNW)
# rm(phosphorusSAG)
# rm(phosphorusUMW)
# rm(phosphorusLowMiss)
#
# colnames(phosphorusAll) <- c("MRB_ID", 'Total_Contributing_Area', 'Total_Upstream_Area', 'Total_Yield', 'Incremental_Load', 'Total_Load', 'HUC8', 'Reach_Name')
#
# # No NA's in nitrogenAll
# anyNA(phosphorusAll)
#
# catch3 <- shapefile('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/2002_sparrow_MRB23457.shp')
#
# # Catches with no corresponding shapefile also have no HUC8 code. There are 12 of them.
# catch3.test <- catch3
# catch3.test@data$Shape <- 1
# catchesNoShape <- join(phosphorusAll, catch3.test@data, by = "MRB_ID", type = 'left')
# nrow(catchesNoShape[ which(is.na(catchesNoShape$Shape)), ])
# rm(catchesNoShape)
# rm(catch3.test)
#
# # One catch in phosphorusAll has two records
# # Remove the second, smaller, record
# notUnique <- vector()
# MRBS <- phosphorusAll$MRB_ID
#
# for (i in 1:length(MRBS)) {
#   thisMRB <- MRBS[i]
#   matches <- which(MRBS == thisMRB)
#   if(length(matches) > 1) notUnique <- c(notUnique, i)
#
#
# }
#
# phosphorusAll[32654, ]
# phosphorusAll[37851, ]
#
# phosphorusAll <- phosphorusAll[c(1:37850, 37852:nrow(phosphorusAll)), ]
#
#
# # Join catches
# catch3@data <- join(catch3@data, phosphorusAll, by = "MRB_ID", type = 'left')
# # There are no catches without modeled values.
# anyNA(catch3$Total_Yield)
#
#
#
# View(catch3@data)
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# catch3 <- spTransform(catch3, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, catch3)
#
# totalYield3.r <- rasterize(catch3, plainUSRaster, "Total_Yield")
# # incLoad3.r <- rasterize(catch3, plainUSRaster, "Incremental_Yield")
#
# writeRaster(totalYield3.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield3', format ='GTiff', overwrite = TRUE, progress = 'text' )
# # writeRaster(incLoad3.r,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/incLoad3', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(phosphorusAll)
# rm(catch3)
# rm(plainUSRaster)
# rm(totalYield3.r)
#
# phosphorus1 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield1.tif')
# phosphorus2 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield2.tif')
# phosphorus3 <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield3.tif')
# phosphorusTogether <- mosaic(phosphorus1, phosphorus2, fun = max)
# phosphorusTogether <- mosaic(phosphorus2, phosphorus3, fun = max)
# writeRaster(phosphorusTogether,'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield_US', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(phosphorus1)
# rm(phosphorus2)
# rm(phosphorus3)
# rm(phosphorusTogether)
#
#
# nitrogen <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Nitrogen Rasters/totalYield_US.tif')
# phosphorus <- raster('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/Phosphorus Rasters/totalYield_US.tif')
# natureServeCountiesSparrow <- natureServeCounties
# natureServeCountiesSparrow <- spTransform(natureServeCounties, CRS(proj4string(nitrogen)))
#
#
# extractSPARROW <- function(i, nutrient) {
#
# # for (i in 1:length(meanYield)) {
#   thisCounty <- natureServeCountiesSparrow[i, ]
#   if(nutrient == "nitrogen") {
#     r <- nitrogen
#   }
#
#   if(nutrient == "phosphorus") {
#     r <- phosphorus
#   }
#   if (!is.null(intersect(thisCounty, r))) {
#     meanYield <- extract(r, thisCounty, fun = mean, na.rm = TRUE)
#   }
#   if (is.null(intersect(thisCounty, r))) {
#     meanYield <- NA
#   }
#
#   return(meanYield)
# }
# meanNYield <- unlist(lapply(1:nrow(natureServeCountiesSparrow), extractSPARROW, nutrient = 'nitrogen'))
# natureServeCountiesSparrow$meanNitrogenYield <- meanNYield
#
# meanPYield <- unlist(lapply(1:nrow(natureServeCountiesSparrow), extractSPARROW, nutrient = 'phosphorus'))
# natureServeCountiesSparrow$meanPhosphorusYield <- meanPYield
#
# write.csv(natureServeCountiesSparrow@data, 'H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/meanYieldsByCounty.csv')
#
# rm(extractSPARROW)
# rm(phosphorus)
# rm(nitrogen)
# rm(natureServeCountiesSparrow)

#### Railroads ####
# railroads <- shapefile('H:/Global Change Program/GIS/Transportation/Railroads/railrdl010g.shp')
#
# natureServeCountiesRailroads <- spTransform(natureServeCounties, CRS(proj4string(railroads)))
#
# getRailKm <- function(i) {
#   thisCounty <- natureServeCountiesRailroads[i, ]
#   if(gIntersects(thisCounty, railroads)) {
#     intersection <- intersect(railroads, thisCounty)
#     railroadKM <- sum(intersection$KILOMETERS, na.omit = TRUE)
#   }
#
#   if(!gIntersects(thisCounty, railroads)) {
#     railroadKM <- NA
#   }
#   return(railroadKM)
# }
#
# natureServeCountiesRailroads$RailKM <- sapply(1:nrow(natureServeCountiesRailroads), getRailKm)
#
# natureServeCountiesRailroads$RailKM[ which(is.na(natureServeCountiesRailroads$RailKM))] <- 0
#
# natureServeCountiesRailroads$RAILROAD_M <- natureServeCountiesRailroads$RailKM * 1000
#
# natureServeCountiesRailroads$railDensity <- natureServeCountiesRailroads$RAILROAD_M / natureServeCountiesRailroads$Shape_Area
#
# write.csv(natureServeCountiesRailroads@data, 'H:/Global Change Program/GIS/Transportation/Railroads/counties_railroad_km.csv')
#
# rm(railroads)
# rm(getRailKm)
# rm(natureServeCountiesRailroads)


#### Oil and gas wells ####
#
# oilAndGas <- shapefile('H:/Global Change Program/GIS/Energy/Oil & Gas Wells - FracTracker.org/US_OG_022014.shp')
#
# natureServeCountiesOilAndGas <- spTransform(natureServeCounties, CRS(proj4string(oilAndGas)))
# # exclude Maryland, North Carolina, and Texas
#
# oilGasCounties <- over(oilAndGas, natureServeCountiesOilAndGas, returnList = FALSE)
# oilGasCounties <- cbind(oilGasCounties, oilAndGas@data)
#
# countOilGas <- function(FIPS_CODE_LONG) {
#   nWells <- nrow(oilGasCounties[ which(oilGasCounties$FIPS_CODE_LONG == FIPS_CODE_LONG), ])
#   if (substring(FIPS_CODE_LONG, 0, 2) == '24') nWells <- NA
#   if (substring(FIPS_CODE_LONG, 0, 2) == '37') nWells <- NA
#   if (substring(FIPS_CODE_LONG, 0, 2) == '48') nWells <- NA
#
#   return(nWells)
# }
#
# natureServeCountiesOilAndGas$OIL_GAS_WELLS <- sapply(natureServeCountiesOilAndGas$FIPS_CODE_LONG, countOilGas)
#
#
# oilGasDf <- natureServeCountiesOilAndGas@data[ ,c('FIPS_CODE_LONG', "OIL_GAS_WELLS")]
#
# oilGasDf
#
# write.csv(oilGasDf, 'H:/Global Change Program/GIS/Energy/Oil & Gas Wells - FracTracker.org/nscOilAndGasWells.csv')
#
# rm(oilGasDf)
# rm(natureServeCountiesOilAndGas)
# rm(countOilGas)
# rm(oilGasCounties)
# rm(oilAndGas)

#### Crude oil pipelines ####
#
# oilPipelines <- shapefile('H:/Global Change Program/GIS/Energy/CrudeOil_Pipelines_US_EIA/CrudeOil_Pipelines_US_201703.shp')
# natureServeCountiesPipelines <- spTransform(natureServeCounties, CRS(proj4string(oilPipelines)))
#
# getPipelineLength <- function(i) {
#   thisCounty <- natureServeCountiesPipelines[i, ]
#   if(gIntersects(thisCounty, oilPipelines)) {
#     intersection <- intersect(oilPipelines, thisCounty)
#     intersection <- spTransform(intersection, CRS(proj4string(oilPipelines)))
#     pipelineLength <- 0
#     for (j in 1:nrow(intersection@data)) {
#       thisline <- gLength(intersection[j, ])
#       pipelineLength <- pipelineLength + thisline
#     }
#   }
#
#   if(!gIntersects(thisCounty, oilPipelines)) {
#     pipelineLength <- NA
#   }
#   return(pipelineLength)
# }
#
# natureServeCountiesPipelines$PipelineLength <- sapply(1:nrow(natureServeCountiesPipelines), getPipelineLength)
#
# natureServeCountiesPipelines$PipelineLength[ which(is.na(natureServeCountiesPipelines$PipelineLength))] <- 0
#
# natureServeCountiesPipelines$PipelineLength.PerArea <- natureServeCountiesPipelines$PipelineLength / natureServeCountiesPipelines$Shape_Area
#
# write.csv(natureServeCountiesPipelines@data, 'H:/Global Change Program/GIS/Energy/CrudeOil_Pipelines_US_EIA/natureServeCounties_PipelineLength.csv')
#
# rm(oilPipelines)
# rm(getPipelineLength)
# rm(natureServeCountiesPipelines)

#### Shale plays ####

# shalePlays <- shapefile('H:/Global Change Program/GIS/Energy/Shale Plays - US Energy Information Administration/US_ShalePlays_EIA_May2011.shp')
# plot(shalePlays)
#
# plainUSRaster <- raster('H:/Global Change Program/GIS/Climate/WORLDCLIM Ver 1pt4 Rel 3/30 arcsec/Elevation - 30 arcsec/elevation.tif')
# plainUSRaster <- plainUSRaster * 0
#
# shaleShape <- spTransform(shalePlays, CRS(projection(plainUSRaster)))
#
# plainUSRaster <- crop(plainUSRaster, shaleShape)
#
# shaleRaster <- rasterize(shaleShape, plainUSRaster)
#
# shaleRaster[!is.na(shaleRaster)] <- 1
#
# shaleRasterMerged <- merge(shaleRaster, plainUSRaster)
#
# writeRaster(shaleRasterMerged, 'H:/Global Change Program/GIS/Energy/Shale Plays - US Energy Information Administration/shalePlays', format ='GTiff', overwrite = TRUE, progress = 'text' )
#
# rm(shalePlays)
# rm(plainUSRaster)
# rm(shaleShape)
# rm(shaleRaster)
# rm(shaleRasterMerged)
#
#
#
# shalePlaysRaster <- raster('H:/Global Change Program/GIS/Energy/Shale Plays - US Energy Information Administration/shalePlays.tif')
#
#
# natureServeCountiesShale <- spTransform(natureServeCounties, CRS(projection(shalePlaysRaster)))
# getShalePlays <- function(FIPS_CODE_LONG) {
#   thisCounty <- natureServeCountiesShale[ which(natureServeCountiesShale$FIPS_CODE_LONG == FIPS_CODE_LONG), ]
#   percentShale <- extract(shalePlaysRaster, thisCounty, fun = mean, na.omit = TRUE)
#   return(percentShale)
# }
#
# PERCENT_SHALE_PLAY <- sapply(natureServeCountiesShale$FIPS_CODE_LONG, getShalePlays)
# PERCENT_SHALE_PLAY <- unlist(PERCENT_SHALE_PLAY)
# fipses <- names(PERCENT_SHALE_PLAY)
# PERCENT_SHALE_PLAY <- as.data.frame(PERCENT_SHALE_PLAY)
# PERCENT_SHALE_PLAY$FIPS_CODE_LONG <-fipses
#
# natureServeCountiesShale@data <- natureServeCountiesShale@data[ , 1:(ncol(natureServeCountiesShale@data) - 1)]
#
# library(plyr)
# natureServeCountiesShale@data <- join(natureServeCountiesShale@data, PERCENT_SHALE_PLAY, by = 'FIPS_CODE_LONG', type = 'left')
# natureServeCountiesShale@data[is.na(natureServeCountiesShale@data$PERCENT_SHALE_PLAY), ] <- 0
#
# write.csv(natureServeCountiesShale@data, 'H:/Global Change Program/GIS/Energy/Shale Plays - US Energy Information Administration/nscShalePlays.csv')
#
# rm(natureServeCountiesShale)
# rm(getShalePlays)
# rm(PERCENT_SHALE_PLAY)
# rm(shalePlaysRaster)
# rm(fipses)


#### Coal mines ####
# ? useful? compare to MRDS
#
# coalMines <- shapefile('H:/Global Change Program/GIS/Energy/CoalMines_US_EIA/CoalMines_US_2014_r3.shp')
# head(coalMines)
# coalMines$STATE_FIPS <- unlist(lapply(coalMines$mstafips, prefix, len = 2, pre = "0"))
# coalMines$CTY_FIPS <- unlist(lapply(coalMines$mctyfips, prefix, len = 3, pre = '0'))
# coalMines$FIPS_CODE_LONG <- paste0(coalMines$STATE_FIPS, coalMines$CTY_FIPS)
#
# countCoalMines <- function(FIPS_CODE_LONG) {
#   nMines <- nrow(coalMines[which(coalMines$FIPS_CODE_LONG == FIPS_CODE_LONG), ])
#   return(nMines)
# }
#
# nCoalMines <- unlist(lapply(natureServeCounties$FIPS_CODE_LONG, countCoalMines))
# natureServeCountiesCoalMines <- cbind(natureServeCounties@data, nCoalMines)
# write.csv(natureServeCountiesCoalMines, 'H:/Global Change Program/GIS/Energy/CoalMines_US_EIA/nscCoalMines.csv')
#
# rm(natureServeCountiesCoalMines)
# rm(coalMines)
# rm(nCoalMines)
# rm(countCoalMines)
#
#
#
#### Add predictor columns to natureServeCounties ####
# Function to standardize columns with mean 0 and standard deviation 1
# Climate change regressions
# Climate stability
# Road lengths and density
# 2010 census - counts and log
# Venter (Human Footprint) layers
# National land cover database
# Active Mines database
# MRDS mines
# Sriram layers
# Crop suitability
# Human appropriation of NPP
# PADUS proportions
# EDDMapS Invasives
# NLCD Percent developed imperviousness
# SPARROW N and P pollution
# Railroads
# Oil and gas wells
# Oil pipelines
# Shale plays
# Coal mines
#
#
#
# # Function to standardize columns with mean 0 and standard deviation 1
# standardize <- function(Xvector) {
#   for (i in 1:length(Xvector)) {
#     Xvector[i] <- (Xvector[i] - mean(Xvector, na.rm=TRUE)) / sd(Xvector, na.rm=TRUE)
#   }
#   return(Xvector)
# }

# # Climate change regressions
# natureServeCountiesClimate <- read.csv('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/climate_change/regressions/natureServeCounties_withAllRegressions.csv')
# natureServeCountiesClimate$FIPS_CODE_LONG <- apply(as.matrix(natureServeCountiesClimate$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# colsToJoin <- c('FIPS_CODE_LONG', 'precipAllYears', 'precipLastThirty', 'tmeanAllYears', 'tmeanLastThirty', 'precipPropAllYears', 'precipPropLastThirty')
# natureServeCountiesClimate <- natureServeCountiesClimate[ , colsToJoin]
# colnames(natureServeCountiesClimate) <- c('FIPS_CODE_LONG', 'preAllYrs', 'preLast30', 'temAllYrs', 'temLast30', 'prePropAllYrs', 'prePropLast30')
#
# natureServeCountiesClimate$temAllYrsStd <- standardize(as.vector(natureServeCountiesClimate$temAllYrs))
# natureServeCountiesClimate$temLast30Std <- standardize(as.vector(natureServeCountiesClimate$temLast30))
#
#
#
# natureServeCountiesClimate$prePropAllYrsStd <- standardize(as.vector(natureServeCountiesClimate$prePropAllYrs))
# natureServeCountiesClimate$prePropAllYrsStdSqrd <- (natureServeCountiesClimate$prePropAllYrsStd) ^ 2
#
# natureServeCountiesClimate$prePropLast30Std <- standardize(as.vector(natureServeCountiesClimate$prePropLast30))
# natureServeCountiesClimate$prePropLast30StdSqrd <- (natureServeCountiesClimate$prePropLast30Std) ^ 2
#
# natureServeCounties@data <- join(natureServeCounties@data, natureServeCountiesClimate, by = 'FIPS_CODE_LONG', type = 'left')
#
#
# rm(natureServeCountiesClimate)
# rm(colsToJoin)

# # Climate stability
#
# kdeStability <- read.csv('H:/Global Change Program/GIS/Climate/PRISM/2 arcmin/biovars_annual/stability/kde_stabilities.csv')
# kdeStability <- kdeStability[ , c('FIPS_CODE_LONG', 'meanStability')]
# colnames(kdeStability) <- c('FIPS_CODE_LONG', 'kdeMeanStab')
# kdeStability$kdeMeanStab <- as.numeric(kdeStability$kdeMeanStab)
# kdeStability$kdeMeanStabStd <- standardize(as.vector(kdeStability$kdeMeanStab))
# kdeStability$FIPS_CODE_LONG <- apply(as.matrix(kdeStability$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
# natureServeCounties@data <- join(natureServeCounties@data, kdeStability, by = 'FIPS_CODE_LONG', type = 'left')
# rm(kdeStability)
#
# # Road lengths and density
# countyRoadLengths <- read.csv('H:/Global Change Program/GIS/Transportation/Roads - USDA/roadsGRASS/county_road_lengths.csv')
# countyRoadLengths$FIPS_CODE_LONG <- as.character(countyRoadLengths$FIPS_CODE_LONG)
# countyRoadLengths$FIPS_CODE_LONG <- unlist(lapply(countyRoadLengths$FIPS_CODE_LONG, prefix, 5, pre = '0'))
# head(countyRoadLengths)
# countyRoadLengths <- countyRoadLengths[,2:3]
# colnames(countyRoadLengths) <- c('TIGERroadLength', 'FIPS_CODE_LONG')
# countyRoadLengths$TIGERroadLength <- unlist(countyRoadLengths$TIGERroadLength)
# natureServeCounties@data <- join(natureServeCounties@data, countyRoadLengths, by = 'FIPS_CODE_LONG', type = 'left')
# #density defined as length of roads (in meters) per unit area of county (Shape_Area - what are the units?)
#
# natureServeCounties$TIGERroadDensity <- natureServeCounties$TIGERroadLength/natureServeCounties$Shape_Area
# natureServeCounties$TIGERroadDensityStd <- standardize(as.vector(natureServeCounties$TIGERroadDensity))
# rm(countyRoadLengths)
#
# # 2010 census
# censusRevised <- read.csv('H:/Global Change Program/GIS/Population/2010 Census Demographic Profile SF/censusClean.csv')
# censusRevised$FIPS_CODE_LONG <- unlist(sapply(censusRevised$FIPS_CODE_LONG, prefix, len = 5, pre ='0'))
# censusRevised$Final_population <- as.numeric(censusRevised$Final_population)
# censusCols <- c('FIPS_CODE_LONG', 'Final_population')
# censusRevised <- censusRevised[,censusCols]
# colnames(censusRevised) <- c('FIPS_CODE_LONG', 'censusPop')
# rm(censusCols)
#
# natureServeCounties@data <- join(natureServeCounties@data, censusRevised, by = 'FIPS_CODE_LONG', type ='left')
# rm(censusRevised)
#
# Census data incomplete matches for Alaska, Puerto Rico.
# noCensus <- natureServeCounties@data[ which(is.na(natureServeCounties@data$Census_population)), ]
# unique(noCensus$STATE_NAME)
# # :Log population
# natureServeCounties$logCensusPop <- log(natureServeCounties$censusPop)
# natureServeCounties$logCensusPopStd <- standardize(as.vector(natureServeCounties$logCensusPop))

# # Venter (Human Footprint) layers
# venter <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Human Footprint/Dryadv2/nscWithVenterMeans.csv')
#
# venter <- venter[ , c('FIPS_CODE_LONG', 'hfiMeans',  'builtMeans', 'nightLightMeans', "pastureMeans",  "cropsMeans", "popdens1990Means", "popdens2010Means", "roadMeans")]
#
# colnames(venter) <- c('FIPS_CODE_LONG', 'venterHFI', 'venterBuilt', 'venterNtLt', 'venterPasture', 'venterCrop', 'venterPop1990s', 'venterPop2010s', "venterRoad")
#
# venter$venterCropPasture <- venter$venterPasture + venter$venterCrop
#
# venter$venterHFIStd <- standardize(as.vector(venter$venterHFI))
# venter$venterBuiltStd <- standardize(as.vector(venter$venterBuilt))
# venter$venterNtLtStd <- standardize(as.vector(venter$venterNtLt))
# venter$venterPastureStd <- standardize(as.vector(venter$venterPasture))
# venter$venterCropStd <- standardize(as.vector(venter$venterCrop))
# venter$venterCropPastureStd <- standardize(as.vector(venter$venterCropPasture))
# venter$venterPop1990sStd <- standardize(as.vector(venter$venterPop1990s))
# venter$venterPop2010sStd <- standardize(as.vector(venter$venterPop2010s))
# venter$venterRoadStd <- standardize(as.vector(venter$venterRoad))
#
# venter$FIPS_CODE_LONG <- apply(as.matrix(venter$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
#
# natureServeCounties@data <- join(natureServeCounties@data, venter, by = 'FIPS_CODE_LONG', type = 'left')
# rm(venter)
#
# # National land cover database
# nlcd <- read.csv('H:/Global Change Program/GIS/Land use/National Land Cover Database 2011/Tallies/nlcdAllTalliesWithFractions.csv')
# nlcd <- nlcd[ , c('FIPS_CODE_LONG', 'Developed.Fraction', 'Cultivated.Fraction', 'Crops.Fraction', 'Pasture.Fraction', 'Grass.Pasture.Fraction')]
#
# colnames(nlcd) <- c('FIPS_CODE_LONG', 'nlcdDev', 'nlcdCultiv', 'nlcdCrops', 'nlcdPasture', 'nlcdGrassPasture')
#
# nlcd$nlcdDevAg <- nlcd$nlcdDev + nlcd$nlcdCultiv
# nlcd$nlcdDevStd <- standardize(as.vector(nlcd$nlcdDev))
# nlcd$nlcdCultivStd <- standardize(as.vector(nlcd$nlcdCultiv))
# nlcd$nlcdCropsStd <- standardize(as.vector(nlcd$nlcdCrops))
# nlcd$nlcdPastureStd <- standardize(as.vector(nlcd$nlcdPasture))
# nlcd$nlcdGrassPastureStd <- standardize(as.vector(nlcd$nlcdGrassPasture))
# nlcd$nlcdDevAgStd <- standardize(as.vector(nlcd$nlcdDevAg))
#
# nlcd$FIPS_CODE_LONG <- apply(as.matrix(nlcd$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
#
# natureServeCounties@data <- join(natureServeCounties@data, nlcd, by = 'FIPS_CODE_LONG')
#
# rm(nlcd)
#
# # Active Mines database
# activeMines <- read.csv('H:/Global Change Program/GIS/Land use/Mining/Active mines in the US/activeMines.csv')
# activeMines <- activeMines[ , 2:3]
#
# colnames(activeMines) <- c('FIPS_CODE_LONG', 'activeMines')
#
# activeMines$FIPS_CODE_LONG <- apply(as.matrix(activeMines$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
#
# natureServeCounties@data <- join(natureServeCounties@data, activeMines, by = "FIPS_CODE_LONG")
#
# natureServeCounties$activeMinesDens <- natureServeCounties$activeMines/natureServeCounties$Shape_Area
# natureServeCounties$activeMinesDensStd <- standardize(as.vector(natureServeCounties$activeMinesDens))
#
# rm(activeMines)
#
#
# # MRDS mines
# mrdsMines <- read.csv('H:/Global Change Program/GIS/Land use/Mining/mrds/mrdsMinesProducers_Plants.csv')
# mrdsMines <- mrdsMines[ , 2:3]
# colnames(mrdsMines) <- c('FIPS_CODE_LONG', 'mrdsMines')
# mrdsMines$FIPS_CODE_LONG <- apply(as.matrix(mrdsMines$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, mrdsMines, by = 'FIPS_CODE_LONG', type = 'left')
# natureServeCounties$mrdsMinesDens <- natureServeCounties$mrdsMines / natureServeCounties$Shape_Area
# natureServeCounties$mrdsMinesDensStd <- standardize(as.vector(natureServeCounties$mrdsMinesDens))
#
#
# rm(mrdsMines)
#
# # Sriram layers
# cars <- read.csv('H:/Global Change Program/GIS/Sriram/cars.csv')
# soilDeg <- read.csv('H:/Global Change Program/GIS/Sriram/SoilDeg.csv')
# proxRoads <- read.csv('H:/Global Change Program/GIS/Sriram/ProxRoads.csv')
# pigs <- read.csv('H:/Global Change Program/GIS/Sriram/Pigs.csv')
# poultry <- read.csv('H:/Global Change Program/GIS/Sriram/Poultry.csv')
#
# sriram <- cbind(cars, soilDeg, proxRoads, pigs, poultry)
# sriram <- sriram[ , c('FIPS_CODE_LONG', 'meanCars', 'meanSoilDeg', 'meanProxRoads', 'meanPigs', 'meanPoultry')]
# colnames(sriram) <- c('FIPS_CODE_LONG', 'sriCars', 'sriSoilDeg', 'sriProxRoads', 'sriPigs', 'sriPoultry')
# sriram$FIPS_CODE_LONG <- apply(as.matrix(sriram$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
#
# sriram$sriCars <- as.numeric(sriram$sriCars)
# sriram$sriSoilDeg <- as.numeric(sriram$sriSoilDeg)
# sriram$sriProxRoads <- as.numeric(sriram$sriProxRoads)
# sriram$sriPigs <- as.numeric(sriram$sriPigs)
# sriram$sriPoultry <- as.numeric(sriram$sriPoultry)
#
#
# sriram$sriCarsStd <- standardize(as.vector(as.numeric(sriram$sriCars)))
# sriram$sriSoilDegStd <- standardize(as.vector(as.numeric(sriram$sriSoilDeg)))
# sriram$sriProxRoadsStd <- standardize(as.vector(as.numeric(sriram$sriProxRoads)))
# sriram$sriPigsStd <- standardize(as.vector(as.numeric(sriram$sriPigs)))
# sriram$sriPoultryStd <- standardize(as.vector(as.numeric(sriram$sriPoultry)))
#
#
#
#
# natureServeCounties@data <- join(natureServeCounties@data, sriram, by = 'FIPS_CODE_LONG', type = 'left')
#
# rm(sriram)
# rm(cars)
# rm(soilDeg)
# rm(proxRoads)
# rm(pigs)
# rm(poultry)
#
#
# # Crop suitability
# cropSuit <- read.csv('H:/Global Change Program/GIS/Land use/Agricultural suitability/meanCropSuit.csv')
# cropSuit <- cropSuit[ ,2:3]
# colnames(cropSuit) <- c('FIPS_CODE_LONG', 'cropSuit')
# cropSuit$cropSuit <- as.numeric(cropSuit$cropSuit)
# cropSuit$cropSuitStd <- standardize(as.vector(cropSuit$cropSuit))
# cropSuit$FIPS_CODE_LONG <- apply(as.matrix(cropSuit$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, cropSuit, by = "FIPS_CODE_LONG", type = "left")
#
# rm(cropSuit)
#
# # Human appropriation of NPP
# hanpp <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Human Appropriation of NPP - SEDAC/meanHANPP.csv')
# hanpp <- hanpp[ ,2:3]
# colnames(hanpp) <- c('FIPS_CODE_LONG', 'hanpp')
# hanpp$hanpp <- as.numeric(hanpp$hanpp)
# hanpp$hanppStd <- standardize(as.vector(hanpp$hanpp))
# hanpp$FIPS_CODE_LONG <- apply(as.matrix(hanpp$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
#
# natureServeCounties@data <- join(natureServeCounties@data, hanpp, by = 'FIPS_CODE_LONG', type = 'left')
#
# rm(hanpp)
#
#
# # PADUS proportions
# padusProportions <- read.csv('H:/Global Change Program/GIS/Protected Areas/PAD-US version 1.4 Combined Feature Class/padusCounties.csv')
# colnames(padusProportions)
# padusProportions <- padusProportions[ , 2:(ncol(padusProportions))]
# for (i in 2:ncol(padusProportions)) {
#   thisColName <- colnames(padusProportions[i])
#   nwColName <- paste0(thisColName, "Std", sep="")
#   thisColStd <- standardize(as.numeric(padusProportions[,i]))
#   thisColStd <- as.data.frame(thisColStd)
#   colnames(thisColStd) <- nwColName
#   padusProportions <- cbind(padusProportions, thisColStd)
# }
# padusProportions$FIPS_CODE_LONG <- apply(as.matrix(padusProportions$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data,padusProportions, by = "FIPS_CODE_LONG", type = "left")
#
# rm(padusProportions)
#
#
#
# # EDDMapS Invasives
# # Not included as of July 2017 - not sure about permissions
# invasives <- read.csv('H:/Global Change Program/GIS/Invasive species/EDDMapS April 12 2017/nsvWithInvasives.csv')
# invasives <- invasives[ , c('FIPS_CODE_LONG', 'Number.of.Invasive.Plants', 'Number.of.Invasive.Diseases', 'Number.of.Invasive.Insects', 'Number.of.Invasive.Wildlife')]
#
# colnames(invasives) <- c('FIPS_CODE_LONG', 'invPlants', 'invDiseases', 'invInsects', 'invAnimals')
# invasives$invPlantsStd <- standardize(as.vector(invasives$invPlants))
# invasives$invDiseasesStd <- standardize(as.vector(invasives$invDiseases))
# invasives$invInsectsStd <- standardize(as.vector(invasives$invInsects))
# invasives$invAnimalsStd <- standardize(as.vector(invasives$invAnimals))
#
#
#
#
# natureServeCounties@data <- join(natureServeCounties@data, invasives, by = 'FIPS_CODE_LONG', type = 'left')
#
# rm(invasives)
#
# # NLCD Percent developed imperviousness
# imperviousness <- read.csv('H:/Global Change Program/GIS/Land use/NCLD Developed Imperviousness 2011/percentImpervious.csv')
# imperviousness <- imperviousness[ , c('FIPS_CODE_LONG', 'Percent.Impervious')]
# colnames(imperviousness) <- c('FIPS_CODE_LONG', 'percentImperv')
# imperviousness$percentImpervStd <- standardize(as.vector(imperviousness$percentImperv))
# imperviousness$FIPS_CODE_LONG <- apply(as.matrix(imperviousness$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, imperviousness, by = "FIPS_CODE_LONG", type = 'left')
#
# rm(imperviousness)
#
#
# # SPARROW N and P pollution
# sparrow <- read.csv('H:/Global Change Program/GIS/Human Impact & Development/Water quality/EPA SPARROW models/SPARROW 2002 Nitrogen and Phosphorous Pollution/meanYieldsByCounty.csv')
# sparrow <- sparrow[ , c('FIPS_CODE_LONG', 'meanNitrogenYield', 'meanPhosphorusYield')]
# colnames(sparrow) <- c('FIPS_CODE_LONG', 'sparrowN', 'sparrowP')
# sparrow$sparrowNP <- sparrow$sparrowN + sparrow$sparrowP
# sparrow$sparrowNStd <- standardize(as.vector(sparrow$sparrowN))
# sparrow$sparrowPStd <- standardize(as.vector(sparrow$sparrowP))
# sparrow$sparrowNPStd <- standardize(as.vector(sparrow$sparrowNP))
# sparrow$FIPS_CODE_LONG <- apply(as.matrix(sparrow$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, sparrow, by = "FIPS_CODE_LONG", type = 'left')
#
# rm(sparrow)
#
# # Railroads
# railroads <- read.csv('H:/Global Change Program/GIS/Transportation/Railroads/counties_railroad_km.csv')
#
#
# railroads <- railroads[ , c('FIPS_CODE_LONG', 'railDensity')]
#
# railroads$railDensityStd <- standardize(railroads$railDensity)
# railroads$FIPS_CODE_LONG <- apply(as.matrix(railroads$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, railroads, by = 'FIPS_CODE_LONG', type ='left')
# rm(railroads)
#
#
# # Oil and gas wells
# oilGas <- read.csv('H:/Global Change Program/GIS/Energy/Oil & Gas Wells - FracTracker.org/nscOilAndGasWells.csv')
# oilGas <- oilGas[, 2:3]
# colnames(oilGas) <- c('FIPS_CODE_LONG', 'nOilGasWells')
# oilGas$nOilGasWellsStd <- standardize(as.vector(oilGas$nOilGasWells))
# oilGas$FIPS_CODE_LONG <- apply(as.matrix(oilGas$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, oilGas, by = 'FIPS_CODE_LONG', type = 'left')
# rm(oilGas)
#
# # Oil pipelines
# oilPipelines <- read.csv('H:/Global Change Program/GIS/Energy/CrudeOil_Pipelines_US_EIA/natureServeCounties_PipelineLength.csv')
# oilPipelines <- oilPipelines[ , c('FIPS_CODE_LONG', 'PipelineLength.PerArea')]
# colnames(oilPipelines) <- c('FIPS_CODE_LONG', 'pipelineLengthPerArea')
# oilPipelines$pipelineLengthPerAreaStd <- standardize(oilPipelines$pipelineLengthPerArea)
# oilPipelines$FIPS_CODE_LONG <- apply(as.matrix(oilPipelines$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, oilPipelines, by = 'FIPS_CODE_LONG', type = 'left')
#
# rm(oilPipelines)
#
# # Shale plays
# shale <- read.csv('H:/Global Change Program/GIS/Energy/Shale Plays - US Energy Information Administration/nscShalePlays.csv')
# shale <- shale[ , c('FIPS_CODE_LONG', 'PERCENT_SHALE_PLAY')]
# colnames(shale) <- c('FIPS_CODE_LONG', 'percentShalePlay')
# shale$percentShalePlayStd <- standardize(as.vector(shale$percentShalePlay))
# shale$FIPS_CODE_LONG <- apply(as.matrix(shale$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, shale, by = 'FIPS_CODE_LONG', type = "left")
# rm(shale)
#
#
# # Coal mines
# coalMines <- read.csv('H:/Global Change Program/GIS/Energy/CoalMines_US_EIA/nscCoalMines.csv')
# coalMines <- coalMines[ , c('FIPS_CODE_LONG', 'nCoalMines')]
#
# coalMines$nCoalMinesStd <- standardize(as.vector(coalMines$nCoalMines))
# coalMines$FIPS_CODE_LONG <- apply(as.matrix(coalMines$FIPS_CODE_LONG), 1, prefix, len = 5, pre = "0")
#
# natureServeCounties@data <- join(natureServeCounties@data, coalMines, by = 'FIPS_CODE_LONG', type = 'left')
#
# rm(coalMines)
#
# rm(standardize)
#
#
#
#### Prepare data for stan
# # Load threat and species occurrence data
#
# speciesThreats <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/Working/56a Master Copy, Assessed Species 2015-08-16 1=past OR present OR future 0=no threat Re-added 7.3 Other modif.csv')
# natureServeCountyOccurrence <- read.csv('//mbgcl02fs/usersdatavol2/CCSD/shared/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/Working/00_NS_mv_CTY_G12ESA_list_201403 - County Occurrences.csv')
#
# # Format FIPS codes correctly
# natureServeCountyOccurrence$FIPS_CODE_LONG <- apply(as.matrix(natureServeCountyOccurrence$FIPS_CD), 1, prefix, len = 5, pre = "0")
#
# # Filter to plants
# speciesThreats <- speciesThreats[ which(speciesThreats$KINGDOM == "Plantae"), ]
# natureServeCountyOccurrence <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$INFORMAL_TAX == "Ferns and relatives" | natureServeCountyOccurrence$INFORMAL_TAX == "Conifers and relatives" | natureServeCountyOccurrence$INFORMAL_TAX == "Hornworts" | natureServeCountyOccurrence$INFORMAL_TAX == "Liverworts" | natureServeCountyOccurrence$INFORMAL_TAX == "Mosses" | natureServeCountyOccurrence$INFORMAL_TAX == "Flowering Plants"), ]
#
#
# # Tally number of species found in each county
# findCountyOccurrences <- function(FIPS_CODE){
#   countyOccurrences <- grep(FIPS_CODE, natureServeCountyOccurrence$FIPS_CODE_LONG)
#   return(as.numeric(length(countyOccurrences)))
# }
#
# natureServeCounties$speciesInCounty <- apply(as.matrix(natureServeCounties$FIPS_CODE_LONG), 1, findCountyOccurrences)
#
#
# # Add conglomerate threats to speciesThreats
# combineThreats <- function(i, threats) {
#   union <- max(speciesThreats[i, threats])
#   return(union)
# }
#
# speciesThreats$c1p0x1x2x3_allDev <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c1p0_residComm', 'c1p1_resident', 'c1p2_comm', 'c1p3_tourDev')))
# speciesThreats$c2p0x1x2x3_allAg <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c2p0_ag', 'c2pt1_crops', 'c2pt2_plantation', 'c2pt3_livestock')))
# speciesThreats$c2p1x3_cropsLive <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c2pt1_crops', 'c2pt3_livestock')))
# speciesThreats$c6p1p0x1x2_allRec <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c6p1p0_rec', 'c6p1p1_hiking', 'c6p1p2_orv')))
# speciesThreats$c6p0x1_allHumanInt <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c6p0_intrusion', 'c6p1p1_hiking', 'c6p1p2_orv', 'c6p2_war', 'c6p3_work')))
# speciesThreats$c9p0x3_pollUnspAg <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c9p0_pollution', 'c9p3_agPollution')))
# speciesThreats$anyThreat <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c(40:90)))
# speciesThreats$allRoads <- unlist(sapply(1:nrow(speciesThreats), combineThreats, threats = c('c4p0_transport', 'c4p1p0_roads', 'c4p1p2_roadConst')))
#
#
# # Option to tally number of species (regardless of endemism) affected by each threat in each county
# countSpeciesThreatened <- function(CTY_FIPS, threat){
#   natureServeCountyIndex <- as.numeric(match(CTY_FIPS, natureServeCounties$FIPS_CODE_LONG))
#   if (natureServeCounties$speciesInCounty[natureServeCountyIndex] == 0) return(0)
#
#   speciesInCounty <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$FIPS_CODE_LONG == CTY_FIPS), ]
#   speciesInCounty <- speciesInCounty$ELEMENT_GLOBAL_ID
#
#   speciesInCountyIndices <- match(speciesInCounty, speciesThreats$ELEMENT_GLOBAL_ID)
#   speciesInCountyIndices <- speciesInCountyIndices[!is.na(speciesInCountyIndices)]
#
#   speciesInCounty <- speciesThreats[speciesInCountyIndices, ]
#
#   speciesThreatsColIndex <- as.numeric(match(threat, colnames(speciesInCounty)))
#
#   speciesThreatenedInCounty <- speciesInCounty[ which(speciesInCounty[speciesThreatsColIndex] == 1), ]
#
#   return(as.numeric(length(speciesThreatenedInCounty$ELEMENT_GLOBAL_ID)))
# }
#
#
# # Tally number of species threatened by a threat
# countSpeciesThreatenedByThreat <- function(threatName) {
#
#   threat <- threatName
#
#   threatCol <- apply(as.matrix(natureServeCounties$FIPS_CODE_LONG), 1, countSpeciesThreatened, threat = threat)
#
#   threatCol <- as.matrix(threatCol)
#
#   colnames(threatCol) <- paste("spThreatenedBy_", threat, sep = "")
#   return(threatCol)
# }
# #
# ## add species affected by each threat in each county to nSC
# # natureServeCountiesT <- natureServeCounties
# # natureServeCountiesT$Index <- row.names(natureServeCountiesT@data)
# # T <- as.data.frame(countSpeciesThreatenedByThreat('c11p3x4_tempPrecip'))
# # T$Index <- row.names(T)
# # natureServeCountiesT@data <- join(natureServeCountiesT@data, T, by = 'Index')
# #
# # natureServeCounties@data <- cbind(natureServeCounties@data, countSpeciesThreatenedByThreat("c4p1p0_roads"))
# #
# # natureServeCounties@data <- cbind(natureServeCounties@data, countSpeciesThreatenedByThreat("c2pt1_crops"))
# #
# # natureServeCounties@data <- cbind(natureServeCounties@data, countSpeciesThreatenedByThreat("c2pt3_livestock"))
#
# # Alternatively, just make a dataframe of nsc + species threatened by all the different threats
#
# threatsList <- colnames(speciesThreats)[c(40:90, 94:length(colnames(speciesThreats)))]
#
# natureServeCountiesThreats <- natureServeCounties@data
# natureServeCountiesThreats$Index <- row.names(natureServeCountiesThreats)
#
# for (i in 1:length(threatsList)) {
#   threatened <- as.data.frame(countSpeciesThreatenedByThreat(threatsList[i]))
#  threatened$Index <- row.names(threatened)
#   natureServeCountiesThreats <- join(natureServeCountiesThreats, threatened, by = 'Index')
#
# }
#
# write.csv(natureServeCountiesThreats, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/natureServeCountiesThreats.csv')
# rm(natureServeCountiesThreats)
# rm(threatsList)
# rm(combineThreats)
# rm(countSpeciesThreatened)
# rm(countSpeciesThreatenedByThreat)
# rm(i)
# rm(findCountyOccurrences)
# rm(natureServeCountyOccurrence)
# rm(speciesThreats)
#
# # Join # of species threatened by each threat with threat layer data
# natureServeCountiesThreats <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/natureServeCountiesThreats.csv')
# colnames(natureServeCountiesThreats)
# natureServeCountiesThreats <- natureServeCountiesThreats[c(17, 19:78)]
# 
# natureServeCounties@data <- join(natureServeCounties@data, natureServeCountiesThreats, by = 'FIPS_CODE_LONG', type = 'left')

natureServeCounties$X_zeros <- 0
natureServeCounties$X_zeros <- as.numeric(natureServeCounties$X_zeros)

saveCols <- colnames(natureServeCounties@data)
write.table(saveCols, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713/cols.txt')

writeOGR(natureServeCounties, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713', layer = 'natureServeCounties', driver = 'ESRI Shapefile', layer_options = 'RESIZE = YES', check_exists = TRUE, overwrite_layer = TRUE)

rm(natureServeCounties)

write.csv(speciesThreats, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713/speciesThreats.csv')


#### CAR stan cascade ####
#### cascade start ####

getFileCheck <- function(threat, predictors, states, ncores, CAR, scramble) {

  predictorString <- predictors[1]
  if (length(predictors) > 1) {
    for (i in 2:length(predictors)) {
      predictorString <- paste(predictorString, predictors[i], sep = "_")
    }
  }
  workPath <- paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master/CAR Analyses/', threat, "_", predictorString, "_st", states, '_scr', scramble, '_CAR', CAR, "/X.tiff")

  return(workPath)
}




setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017')
source('Stan function.R')

#8 gb
runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

# 32gb 3
runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32 gb 3
runStanModel('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32gb4
runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)


#### cascade 2 ####

#8 gb
file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

# 32gb 3
file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32 gb 3
file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32gb4
file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)



while(!(file.exists(file))) {

  # just wait

}

setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017')
source('Stan function.R')


#8 gb
#file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

# 32gb 3
#file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32 gb 3
# file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
runStanModel('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32gb4
#file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)
runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)



#### cascade 3 ####
#8 gb
#file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
# runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

# 32gb 3
#file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32 gb 3
# file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('allRoads', c('TIGERroadDensityStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32gb4
#file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)
file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
#runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)


while(!(file.exists(file))) {

  # just wait

}

setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017')
source('Stan function.R')

#8 gb
#file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#file <- getFileCheck(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
runStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

# 32gb 3
#file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#file <- getFileCheck(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
runStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32 gb 3
# file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE)
#file <- getFileCheck('allRoads', c('TIGERroadDensityStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
runStanModel('allRoads', c('TIGERroadDensityStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)

#32gb4
#file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)
#file <- getFileCheck('c11p3x4_tempPrecip', c('temAllYrsStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE)
runStanModel('c11p3x4_tempPrecip', c('temAllYrsStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE)


#
# scrambledNSC <- natureServeCounties
# scrambledNSC$index <- c(1:nrow(scrambledNSC@data))
# dat <- scrambledNSC@data[, c(1:3, 5:15, 17:ncol(scrambledNSC@data))]
# dat$FIPS_CODE_LONG_OLD <- natureServeCounties$FIPS_CODE_LONG
# fips <- scrambledNSC@data[,c('FIPS_CODE_LONG','STATE_NAME', 'index')]
# set.seed(10)
# fips$index <- sample.int(nrow(fips), nrow(fips), replace = FALSE)
# fips <- join(fips, dat, by = 'index', type='left')
# scrambledNSC@data <- fips
#
#
# key <- fips[,c('FIPS_CODE_LONG', 'FIPS_CODE_LONG_OLD')]
#
# scrambledOcc <- natureServeCountyOccurrence
# colnames(scrambledOcc) <- c(colnames(scrambledOcc)[1:18], 'FIPS_CODE_LONG_OLD')
# scrambledOcc <- join(scrambledOcc, key, by = 'FIPS_CODE_LONG_OLD', type = 'left')
#
# # runScrambledCARModel <- function(threat, predictors, states, ncores){
#
#   # get X and list of counties that have non-NA threat values
#   #if (states == TRUE) colsToUse <- c('FIPS_CODE_LONG', predictors, 'STATE_NAME')
#   #if (states == FALSE) colsToUse <- c('FIPS_CODE_LONG', predictors)
#   colsToUse <- c('FIPS_CODE_LONG', predictors, 'STATE_NAME')
#   predictorDf <- scrambledNSC[ , colsToUse]
#   colnames(predictorDf@data) <- c('FIPS_CODE_LONG', 'X','STATE_NAME')
#   for (i in 1:length(predictors)) {
#     thisCol <- predictorDf@data[ , i + 1]
#     predictorDf <- predictorDf[which(!is.na(thisCol)), ]
#   }
#
#   states.num <- unique(predictorDf$STATE_NAME)
#   states.num <- as.data.frame(states.num)
#   states.num$state.hyperP <- c(1:nrow(states.num))
#   colnames(states.num) <- c('STATE_NAME', "STATEFP")
#   predictorDf@data <- join(predictorDf@data, states.num, by = "STATE_NAME")
#
#
#
#   nCounties <- nrow(predictorDf)
#
#   predictorCols <- 2:(length(predictors) + 1)
#
#
#   X <- predictorDf@data[,predictorCols]
#
#
#   X <- as.matrix(X)
#
#   allSpecies <- unique(speciesThreats$ELEMENT_GLOBAL_ID)
#
#   occMatrix <- matrix(data=0, nrow = length(allSpecies), ncol = nrow(predictorDf))
#   for (i in 1:nrow(occMatrix)) {
#     thisSpecies <- allSpecies[i]
#     speciesCounties <- scrambledOcc[ which(scrambledOcc$ELEMENT_GLOBAL_ID == thisSpecies), 'FIPS_CODE_LONG']
#     for (j in 1:length(speciesCounties)){
#       countyIndex <- as.numeric(match(speciesCounties[j], predictorDf$FIPS_CODE_LONG))
#       occMatrix[i,countyIndex] <- 1
#     }
#   }
#
#   speciesZeros <- vector(length=length(allSpecies))
#   for (i in 1:length(allSpecies)) {
#     speciesZeros[i] <- sum(occMatrix[i, ])
#   }
#
#   speciesToKeep <- cbind(allSpecies, speciesZeros)
#   speciesToKeep <- as.data.frame(speciesToKeep)
#   speciesToKeep$indices <- row.names(speciesToKeep)
#   speciesToKeep <- speciesToKeep[ which(speciesToKeep$speciesZeros > 0), ]
#
#   speciesIndices <- as.vector(as.numeric(speciesToKeep$indices))
#   nSpecies <- length(speciesIndices)
#
#   occMatrix <- occMatrix[speciesIndices, ]
#   predictorDf$Richness <- NA
#   for (i in 1:nCounties) {
#     predictorDf$Richness[i] <- sum(occMatrix[,i])
#   }
#   # at this point occMatrix has a row for each species that is not excluded from this analysis and a column for each county
#   # CAR needs natureServeCounties to have a column for each species that says whether it is present in a county
#   SppCols <- vector(length = length(speciesIndices))
#
#   for (i in 1:length(speciesIndices)) {
#     occVect <- as.vector(occMatrix[i, ])
#     predictorDf@data <- cbind(predictorDf@data, occVect)
#     colnames(predictorDf@data) <- c(colnames(predictorDf@data[1:ncol(predictorDf@data) - 1]), paste0('spp', speciesIndices[i]))
#     SppCols[i] <- paste0('spp', speciesIndices[i])
#
#   }
#
#
#
#   threatIndex <- as.numeric(match(threat, colnames(speciesThreats)))
#
#   EndSpec <- speciesToKeep
#   EndSpec$end <- 0
#   for(i in 1:nrow(speciesToKeep)){
#     thisSpec <- speciesThreats[ which(speciesThreats$ELEMENT_GLOBAL_ID == speciesToKeep$allSpecies[i]), ]
#     if (thisSpec[1,threatIndex] == 1) EndSpec$end[i] <- 1
#   }
#
#   EndSpec <- as.matrix(EndSpec$end)
#
#   if (states == FALSE) {
#     nHyperP = 1
#     HyperPAssign = vector( length = nCounties)
#
#
#     for(i in 1:nCounties) {
#
#       HyperPAssign[i] <- 1
#     }
#   }
#
#   if (states == TRUE) {
#     nHyperP = as.numeric(length(unique(predictorDf$STATEFP)))
#     HyperPAssign = predictorDf$STATEFP
#   }
#
#   HyperPAssign <- as.numeric(HyperPAssign)
#
#   setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June')
#
#   source("shapeLopodData.R")
#   library(rgeos)
#   library(raster)
#   library(ggplot2)
#   library(slam)
#
#
#   AdMatrixData = shapeLopodData(predictorDf, fieldN = "Richness", fieldY = "Richness", keepFields = TRUE, Adjacency = TRUE)
#
#
#   #Load Stan
#   library(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())
#
#   #Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
#   #X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)
#
#   #Make sure there are objects nCounties and nSpecies with the number of species and counties respectively
#
#   AllIncludedFeatures = rbind(AdMatrixData$geoInfo$sampledId,AdMatrixData$geoInfo$notSampledId)
#   AllIncludedFeatures = AllIncludedFeatures[order(AllIncludedFeatures$cellStan),]
#
#   #Create x predictors Matrix
#   X = as.matrix(cbind(AdMatrixData$geoDataObject@data[AllIncludedFeatures$featureShape,"X"]))
#
#
#   whichEnd = which(EndSpec == 1)
#   whichNotEnd  = which(EndSpec == 0)
#
#   nEnd = length(whichEnd)
#   nNotEnd = length(whichNotEnd)
#
#   #Field in shapefile with state
#   nHyperP = length(unique(AdMatrixData$geoDataObject@data$STATEFP))
#   HyperPAssign = as.numeric(factor(as.character(AdMatrixData$geoDataObject@data$STATEFP[AllIncludedFeatures$featureShape])))
#
#   spOccMat = t(AdMatrixData$geoDataObject@data[AdMatrixData$geoInfo$sampledId$featureShape,SppCols])
#
#   # create directory and save data
#   predictorString <- predictors[1]
#   if (length(predictors) > 1) {
#     for (i in 2:length(predictors)) {
#       predictorString <- paste(predictorString, predictors[i], sep = "_")
#     }
#   }
#   workPath <- paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled/', threat, "_", predictorString, "_st", states)
#   if (dir.exists(workPath)) setwd(workPath)
#   if (!dir.exists(workPath)) {
#     dir.create(workPath)
#     setwd(workPath)
#   }
#
#   saveRows <- c('nCounties', 'nSpecies', 'nEnd', 'nNotEnd', 'nHyperP', 'K')
#   saveItems <- as.data.frame(saveRows)
#   saveItems$Values <- c(nCounties, nSpecies, nEnd, nNotEnd, nHyperP, dim(X)[2])
#   write.csv(saveItems, "nValues.csv", row.names=FALSE)
#   rm(saveRows)
#   rm(saveItems)
#
#   write.csv(occMatrix, "occMatrix.csv", row.names = FALSE)
#   write.table(whichEnd, "whichEnd.csv", row.names = FALSE)
#   write.table(whichNotEnd, "whichNotEnd.csv", row.names = FALSE)
#   write.table(HyperPAssign, "HyperPAssign.csv", row.names = FALSE)
#   write.csv(X, "X.csv", row.names = FALSE)
#   write.csv(speciesToKeep, "speciesToKeep.csv", row.names=FALSE)
#   write.csv(predictorDf@data, "predictorDf.csv", row.names = FALSE)
#
#   setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June')
#
#   #Names of the list elements should be the same (left side of the =)
#   stanData = list( nCounties = length(AdMatrixData$geoInfo$sampledId$featureShape),
#                    sampledId = AdMatrixData$geoInfo$sampledId$cellStan,
#                    nNotSampled = length(AdMatrixData$geoInfo$notSampledId$featureShape),
#                    notSampledId = AdMatrixData$geoInfo$notSampledId$cellStan,
#                    n = length(AdMatrixData$geoInfo$sampledId$featureShape) + length(AdMatrixData$geoInfo$notSampledId$featureShape) ,
#                    W_n = dim(AdMatrixData$geoInfo$W_sparse)[1] ,
#                    W_sparse = AdMatrixData$geoInfo$W_sparse ,
#                    D_sparse = AdMatrixData$geoInfo$D_sparse ,
#                    lambda = AdMatrixData$geoInfo$lambda_sparse  ,
#                    nSpecies = nEnd+nNotEnd,
#                    nEnd = nEnd,
#                    nNotEnd = nNotEnd,
#                    spOccMat = spOccMat,#Species by county Matrix (1 and 0)
#                    endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
#                    notEndSpp = whichNotEnd, #Which species are NOT endangered (numbers matching spOccMat)
#                    nHyperP = nHyperP, #Number of categories (States)
#                    HyperPAssign = HyperPAssign, #Vector of length counties assigning it to a state (numeric)
#                    K = dim(X)[2],
#                    x_pred = X
#
#   )
#
#
#   #This line loads the compiled C++ Stan Model, if all files are in the same working directory, you shouldnt need to compile again, unless you change the .stan file
#   StanModel = stan_model(file = "countyEndangered_NotEndemics_CAR.stan" )
#
#   #Runs the MCMC model
#   FitModel = sampling(StanModel,
#                       data = stanData,              # named list of data
#                       chains = 4,                   # number of Markov chains
#                       warmup = 1000,               # number of warmup iterations per chain
#                       iter = 2000,                 # total number of iterations per chain (includes warm-up)
#                       cores = ncores,                    # number of cores
#                       refresh = 10                 # show progress every 'refresh' iterations
#   )
#
#   #Calls coeff for predictors in matrix X and intersect
#   summ1 <-summary(FitModel, pars = c("a",paste("b[",1:dim(X)[2],"]", sep="")))$summary
#   if(is.null(summ1)) return()
#   done <- matrix(data = 1, nrow = 1, ncol = 1)
#   write.table(done, paste0(workPath, '/done.csv'))
#   rm(done)
#   #Calls r squareds
#   summ2 <- summary(FitModel, pars = c("r_sq","r_sq_justX","r_sq_notGeo" ))$c_summary[,,2] #[,"50%"]
#
#   save(FitModel, file= paste0(workPath, "/FitModel_NoEndemics.Rdata"))
#   resultsShape = AdMatrixData$geoDataObject
#   resultsShape@data[AllIncludedFeatures$featureShape,"state_coeff"] = summary(FitModel, pars =paste("a_cat[",1:nHyperP,"]", sep=""))$c_summary[,,2][paste("a_cat[",HyperPAssign,"]", sep=""),"50%"]
#   resultsShape@data[AdMatrixData$geoInfo$sampledId$featureShape,"p_est"] = summary(FitModel, pars = paste("p[",1:length(AdMatrixData$geoInfo$sampledId$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc"] = summary(FitModel, pars = paste("calc_p[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_justX"] = summary(FitModel, pars = paste("calc_p_justX[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_notGeo"] = summary(FitModel, pars = paste("calc_p_notGeo[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"geo_effect"] = summary(FitModel, pars = paste("geo_effect[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data <- resultsShape@data[ , c('FIPS_CODE_LONG', 'X', 'STATE_NAME', 'STATEFP', 'Richness', 'detections', 'FeatureID', 'state_coeff', 'p_est', 'p_calc', 'p_calc_justX', 'p_calc_notGeo' ,'geo_effect')]
#
#   blues <- brewer.pal(9, 'Blues')
#   blues <- c(blues, '#000000')
#   blues <-  cbind(blues, c(0:9))
#   colnames(blues) <- c('Colour', 'pColour')
#   blues <- as.data.frame(blues)
#
#   blues_p_est <- blues
#   blues_p_calc<- blues
#   blues_p_calc_notGeo <- blues
#   blues_p_calc_justX <- blues
#
#   colnames(blues_p_est) <- c('p_estColour', 'p_estVal')
#   colnames(blues_p_calc) <- c('p_calcColour', 'p_calcVal')
#   colnames(blues_p_calc_notGeo) <- c('p_calc_notGeoColour', 'p_calc_notGeoVal')
#   colnames(blues_p_calc_justX) <- c('p_calc_justXColour', 'p_calc_justXVal')
#
#   resultsShape$p_estVal <- floor(as.numeric(resultsShape$p_est) /0.1)
#   resultsShape$p_calcVal <- floor(as.numeric(resultsShape$p_calc) / 0.1)
#   resultsShape$p_calc_notGeoVal <- floor(as.numeric(resultsShape$p_calc_notGeo) / 0.1)
#   resultsShape$p_calc_justXVal <- floor(as.numeric(resultsShape$p_calc_justX) / 0.1)
#
#
#   resultsShape@data <- join(resultsShape@data, blues_p_est, by = 'p_estVal', type = 'left')
#   resultsShape@data <- join(resultsShape@data, blues_p_calc, by = 'p_calcVal', type = 'left')
#   resultsShape@data <- join(resultsShape@data, blues_p_calc_notGeo, by = 'p_calc_notGeoVal', type = 'left')
#   resultsShape@data <- join(resultsShape@data, blues_p_calc_justX, by = 'p_calc_justXVal', type = 'left')
#
#
#   png('p_est.png')
#   plot(resultsShape, col = resultsShape$p_estColour, border = NA)
#   dev.off()
#   png('p_calc.png')
#   plot(resultsShape, col = resultsShape$p_calcColour, border = NA)
#   dev.off()
#   png('p_calc_notGeo.png')
#   plot(resultsShape, col = resultsShape$p_calc_notGeoColour, border = NA)
#   dev.off()
#   png('p_calc_justX.png')
#   plot(resultsShape, col = resultsShape$p_calc_justXColour, border = NA)
#   dev.off()
#
#    shapefile(resultsShape, filename = paste0(workPath, "/countyOutput_CAR.shp"), overwrite = TRUE)
#
#
#   results <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/results.csv')
#   results <- results[,2:29]
#   resultnames <- c('threat', 'b1', 'b2', 'states', 'scrambled', 'a.rhat', 'a.mean', 'a.2pt5', 'a.97pt5', 'b1.rhat', 'b1.mean', 'b1.2pt5', 'b1.95pt5', 'b1.sig', 'b2.rhat', 'b2.mean', 'b2.2pt5', 'b2.97pt5', 'b2.sig', 'rsquared.mean', 'rsquared.2pt5', 'rsquared.97.5', 'rsquared.justX.mean', 'rsquared.justX.2pt5', 'rsquared.justX.97pt5', 'rsquared.notGeo.mean', 'rsquared.notGeo.2pt5', 'rsquared.notGeo.97.5')
#
#   threat <- threat
#   states <- states
#   scrambled <- TRUE
#   b1 <- predictors[1]
#   b2 <- NA
#   if (length(predictors) > 1) b2 <- predictors[2]
#
#
#   a.rhat <- summ1[1,10]
#   a.mean <- summ1[1,1]
#   a.2pt5<- summ1[1,4]
#   a.97pt5 <- summ1[1,8]
#   a.sig <- (a.2pt5 * a.97pt5 > 0)
#
#   b1.rhat <- summ1[2,10]
#   b1.mean <- summ1[2,1]
#   b1.2pt5 <- summ1[2,4]
#   b1.97pt5 <- summ1[2,8]
#   b1.sig <- (b1.2pt5 * b1.97pt5 > 0)
#
#   b2.rhat <- NA
#   b2.mean <- NA
#   b2.2pt5 <-NA
#   b2.97pt5 <- NA
#   b2.sig <- NA
#
#   if (!is.na(b2)) {
#     b2.rhat <- summ1[3,10]
#     b2.mean <- summ1[3,1]
#     b2.2pt5 <- summ1[3,4]
#     b2.97pt5 <- summ1[3,8]
#     b2.sig <- (b2.2pt5 * b2.97pt5 > 0)
#   }
#
#
#   rsquared.mean <- summ2[1,1]
#   rsquared.2pt5 <- summ2[1,3]
#   rsquared.97pt5 <- summ2[1,7]
#
#
#   rsquared.justX.mean <- summ2[2,1]
#   rsquared.justX.2pt5 <- summ2[2,3]
#   rsquared.justX.97pt5 <- summ2[2,7]
#
#
#   rsquared.notGeo.mean <- summ2[3,1]
#   rsquared.notGeo.2pt5 <- summ2[3,3]
#   rsquared.notGeo.97pt5 <- summ2[3,7]
#
#   thisrow <- c(threat, b1, b2, states, scrambled, a.rhat, a.mean, a.2pt5, a.97pt5, b1.rhat, b1.mean, b1.2pt5, b1.97pt5, b1.sig, b2.rhat, b2.mean, b2.2pt5, b2.97pt5, b2.sig, rsquared.mean, rsquared.2pt5, rsquared.97pt5, rsquared.justX.mean, rsquared.justX.2pt5, rsquared.justX.97pt5, rsquared.notGeo.mean, rsquared.notGeo.2pt5, rsquared.notGeo.97pt5)
#   if (length(thisrow) < 28) return()
#   thisrow <- matrix(data = thisrow, nrow = 1, ncol = 28, dimnames = list('NA', resultnames))
#   results <- rbind(results, thisrow)
#   write.csv(results, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/results.csv')
#
#   #shapefile(resultsShape, filename = paste0(workPath, "/countyOutput_CAR.shp"))
#
# }
#
#
# runALot <- function(threat, predictors, states, ncores, scramble) {
#   predictorString <- predictors[1]
#   if (length(predictors) > 1) {
#     for (i in 2:length(predictors)) {
#       predictorString <- paste(predictorString, predictors[i], sep = "_")
#     }
#   }
#   if (scramble) workPath <- paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled/', threat, "_", predictorString, "_st", states)
#   if(!scramble) workPath <-   paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/', threat, "_", predictorString, "_st", states)
#
#   if (dir.exists(workPath)) setwd(workPath)
#   if (!dir.exists(workPath)) {
#     dir.create(workPath)
#     setwd(workPath)
#   }
#
#   done <- file.exists(paste0(workPath, '/done.csv'))
#   while(!done) {
#     if (scramble) runScrambledCARModel(threat, predictors, states, ncores)
#     if (!scramble) runCARModel(threat, predictors, states, ncores)
#     done <- file.exists(paste0(workPath, '/done.csv'))
#   }
#
# }
# #


#### ANOVAS ####
#
#
# testStates <- function(predictor){
#
#   thisDf <- as.data.frame(natureServeCounties[, c('STATE_NAME', predictor)])
#
#   thisDf <- na.omit(thisDf)
#
#   colnames(thisDf) <- c('STATE_NAME', 'X')
#
#   this.aov <- aov(X ~ STATE_NAME, data = thisDf)
#   this.aov.summ <- summary(this.aov)
#   this.aov.sig <- as.numeric(unlist(this.aov.summ)[9])
#
#
#   this.lm <- lm(X ~ STATE_NAME, data = thisDf)
#   this.lm.summ <- summary(this.lm)
#   this.lm.mr2 <- this.lm.summ[[8]]
#   this.lm.adjr2 <- this.lm.summ[[9]]
#
#   stats <- paste(predictor, this.aov.sig, this.lm.mr2, this.lm.adjr2, "\n", sep = ",")
#   cat(stats, file = 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/States_predictors_comparisons.csv', append = TRUE)
#
# }
#
# # statshead <- paste('Predictor', 'ANOVA P(>F)', 'lm Multiple R-squared', 'lm Adjusted R-squared', sep = ",")
# # writeLines(statshead, con = 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/States_predictors_comparisons.csv', sep = "\n")
#
# predictors <- colnames(natureServeCounties@data)[17:104]
#
# for (i in 1:length(predictors)) {
#   testStates(predictors[i])
# }
#
# stateComparisons <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/States_predictors_comparisons.csv', row.names = NULL)
# tops <- colnames(stateComparisons)
# tops <- tops[2:5]
# stateComparisons <- stateComparisons[,1:4]
# colnames(stateComparisons) <- tops
# rm(tops)
#
# stateComparisons$ANOVA.Sig <- (stateComparisons$ANOVA.P..F. <= 0.05)
# # Old CAR code
# readCARModel('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c2pt1_crops_cropSuitStd_stTRUE')
#
# #save(FitModel, file= paste("./",format(Sys.time(), "%b%d%Y"),"FitModel_NoEndemics.RData",sep =""))
#
# load('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c1p0x1x2x3_allDev_nlcdDevStd_stTRUE/FitModel_NoEndemics.Rdata')
#
# #### Save shapefles with results data - for models run before June 22 2017 ####
#
# saveShape <- function(threat, predictors, states, ncores){
#
#   # get X and list of counties that have non-NA threat values
#   #if (states == TRUE) colsToUse <- c('FIPS_CODE_LONG', predictors, 'STATE_NAME')
#   #if (states == FALSE) colsToUse <- c('FIPS_CODE_LONG', predictors)
#   colsToUse <- c('FIPS_CODE_LONG', predictors, 'STATE_NAME')
#   predictorDf <- natureServeCounties[ , colsToUse]
#   colnames(predictorDf@data) <- c('FIPS_CODE_LONG', 'X','STATE_NAME')
#   for (i in 1:length(predictors)) {
#     thisCol <- predictorDf@data[ , i + 1]
#     predictorDf <- predictorDf[which(!is.na(thisCol)), ]
#   }
#
#   states.num <- unique(predictorDf$STATE_NAME)
#   states.num <- as.data.frame(states.num)
#   states.num$state.hyperP <- c(1:nrow(states.num))
#   colnames(states.num) <- c('STATE_NAME', "STATEFP")
#   predictorDf@data <- join(predictorDf@data, states.num, by = "STATE_NAME")
#
#
#
#   nCounties <- nrow(predictorDf)
#
#   predictorCols <- 2:(length(predictors) + 1)
#
#
#   X <- predictorDf@data[,predictorCols]
#
#
#   X <- as.matrix(X)
#
#   allSpecies <- unique(speciesThreats$ELEMENT_GLOBAL_ID)
#
#   occMatrix <- matrix(data=0, nrow = length(allSpecies), ncol = nrow(predictorDf))
#   for (i in 1:nrow(occMatrix)) {
#     thisSpecies <- allSpecies[i]
#     speciesCounties <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$ELEMENT_GLOBAL_ID == thisSpecies), 'FIPS_CODE_LONG']
#     for (j in 1:length(speciesCounties)){
#       countyIndex <- as.numeric(match(speciesCounties[j], predictorDf$FIPS_CODE_LONG))
#       occMatrix[i,countyIndex] <- 1
#     }
#   }
#
#   speciesZeros <- vector(length=length(allSpecies))
#   for (i in 1:length(allSpecies)) {
#     speciesZeros[i] <- sum(occMatrix[i, ])
#   }
#
#   speciesToKeep <- cbind(allSpecies, speciesZeros)
#   speciesToKeep <- as.data.frame(speciesToKeep)
#   speciesToKeep$indices <- row.names(speciesToKeep)
#   speciesToKeep <- speciesToKeep[ which(speciesToKeep$speciesZeros > 0), ]
#
#   speciesIndices <- as.vector(as.numeric(speciesToKeep$indices))
#   nSpecies <- length(speciesIndices)
#
#   occMatrix <- occMatrix[speciesIndices, ]
#   predictorDf$Richness <- NA
#   for (i in 1:nCounties) {
#     predictorDf$Richness[i] <- sum(occMatrix[,i])
#   }
#   # at this point occMatrix has a row for each species that is not excluded from this analysis and a column for each county
#   # CAR needs natureServeCounties to have a column for each species that says whether it is present in a county
#   SppCols <- vector(length = length(speciesIndices))
#
#   for (i in 1:length(speciesIndices)) {
#     occVect <- as.vector(occMatrix[i, ])
#     predictorDf@data <- cbind(predictorDf@data, occVect)
#     colnames(predictorDf@data) <- c(colnames(predictorDf@data[1:ncol(predictorDf@data) - 1]), paste0('spp', speciesIndices[i]))
#     SppCols[i] <- paste0('spp', speciesIndices[i])
#
#   }
#
#
#
#   threatIndex <- as.numeric(match(threat, colnames(speciesThreats)))
#
#   EndSpec <- speciesToKeep
#   EndSpec$end <- 0
#   for(i in 1:nrow(speciesToKeep)){
#     thisSpec <- speciesThreats[ which(speciesThreats$ELEMENT_GLOBAL_ID == speciesToKeep$allSpecies[i]), ]
#     if (thisSpec[1,threatIndex] == 1) EndSpec$end[i] <- 1
#   }
#
#   EndSpec <- as.matrix(EndSpec$end)
#
#   if (states == FALSE) {
#     nHyperP = 1
#     HyperPAssign = vector( length = nCounties)
#
#
#     for(i in 1:nCounties) {
#
#       HyperPAssign[i] <- 1
#     }
#   }
#
#   if (states == TRUE) {
#     nHyperP = as.numeric(length(unique(predictorDf$STATEFP)))
#     HyperPAssign = predictorDf$STATEFP
#   }
#
#   HyperPAssign <- as.numeric(HyperPAssign)
#
#   setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June')
#
#   source("shapeLopodData.R")
#   library(rgeos)
#   library(raster)
#   library(ggplot2)
#   library(slam)
#
#
#   AdMatrixData = shapeLopodData(predictorDf, fieldN = "Richness", fieldY = "Richness", keepFields = TRUE, Adjacency = TRUE)
#
#
#   #Load Stan
#   library(rstan)
#   rstan_options(auto_write = TRUE)
#   options(mc.cores = parallel::detectCores())
#
#   #Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
#   #X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)
#
#   #Make sure there are objects nCounties and nSpecies with the number of species and counties respectively
#
#   AllIncludedFeatures = rbind(AdMatrixData$geoInfo$sampledId,AdMatrixData$geoInfo$notSampledId)
#   AllIncludedFeatures = AllIncludedFeatures[order(AllIncludedFeatures$cellStan),]
#
#   #Create x predictors Matrix
#   X = as.matrix(cbind(AdMatrixData$geoDataObject@data[AllIncludedFeatures$featureShape,"X"]))
#
#
#   whichEnd = which(EndSpec == 1)
#   whichNotEnd  = which(EndSpec == 0)
#
#   nEnd = length(whichEnd)
#   nNotEnd = length(whichNotEnd)
#
#   #Field in shapefile with state
#   nHyperP = length(unique(AdMatrixData$geoDataObject@data$STATEFP))
#   HyperPAssign = as.numeric(factor(as.character(AdMatrixData$geoDataObject@data$STATEFP[AllIncludedFeatures$featureShape])))
#
#   spOccMat = t(AdMatrixData$geoDataObject@data[AdMatrixData$geoInfo$sampledId$featureShape,SppCols])
#
#   # create directory and save data
#   predictorString <- predictors[1]
#   if (length(predictors) > 1) {
#     for (i in 2:length(predictors)) {
#       predictorString <- paste(predictorString, predictors[i], sep = "_")
#     }
#   }
#   workPath <- paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/', threat, "_", predictorString, "_st", states)
#   if (dir.exists(workPath)) setwd(workPath)
#   if (!dir.exists(workPath)) {
#     dir.create(workPath)
#     setwd(workPath)
#   }
#
#   load(paste0(workPath, "/FitModel_NoEndemics.Rdata"))
#
#
#   resultsShape = AdMatrixData$geoDataObject
#   resultsShape@data[AllIncludedFeatures$featureShape,"state_coeff"] = summary(FitModel, pars =paste("a_cat[",1:nHyperP,"]", sep=""))$c_summary[,,2][paste("a_cat[",HyperPAssign,"]", sep=""),"50%"]
#   resultsShape@data[AdMatrixData$geoInfo$sampledId$featureShape,"p_est"] = summary(FitModel, pars = paste("p[",1:length(AdMatrixData$geoInfo$sampledId$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc"] = summary(FitModel, pars = paste("calc_p[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_justX"] = summary(FitModel, pars = paste("calc_p_justX[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_notGeo"] = summary(FitModel, pars = paste("calc_p_notGeo[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data[AllIncludedFeatures$featureShape,"geo_effect"] = summary(FitModel, pars = paste("geo_effect[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,2][,"50%"]
#   resultsShape@data <- resultsShape@data[ , c(1:5, 2186:2194)]
#
#   shapefile(resultsShape, "countyOutput_CAR.shp", overwrite = TRUE)
# }

# library(rstan)
# readCARModel <- function(workPath) {
#   results <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/results.csv')
#   results <- results[,2:29]
#   scrambled <- FALSE
#   split <- strsplit(workPath, "Analyses/")
#   split <- unlist(split)
#   if (substr(split[2],1,2) == 'sc') {
#     split <- strsplit(split[2], 'scrambled/')
#     scrambled <- TRUE
#   }
#   split <- unlist(split)[2]
#   split <- strsplit(split, "_")
#   split <- unlist(split)
#
#   if (split[1] == 'anyThreat') {
#     threat <- split[1]
#     if (split[2] =='X') {
#       b1 <- 'X_zeros'
#     } else {
#       b1 <- split[2]
#     }
#
#     if(length(split) > 3) {
#       if(split[2]!='X')
#         b2 <- split[3]
#       states <- split[4]
#     } else {
#       b2 <- NA
#     }
#   } else {
#     threat <- paste0(split[1], '_', split[2])
#     if(split[3] =='X') {
#       b1 <- 'X_zeros'
#     } else {
#       b1 <- split[3]
#     }
#     if (length(split) > 4) {
#       if ( split[3] != 'X') {
#         b2 <- split[4]
#       }
#     } else {
#       b2 <- NA
#     }
#   }
#
#   states <- split[length(split)]
#
#
#   setwd(workPath)
#
#   savedNs <- read.csv("nValues.csv")
#   nCounties <- as.numeric(savedNs[1,2])
#   nSpecies <- as.numeric(savedNs[2,2])
#   nEnd <- as.numeric(savedNs[3,2])
#   nNotEnd <- as.numeric(savedNs[4,2])
#   nHyperP <- as.numeric(savedNs[5,2])
#   K <- as.numeric(savedNs[6,2])
#   rm(savedNs)
#
#   occMatrix <- read.csv("occMatrix.csv")
#   whichEnd <- read.csv("whichEnd.csv")
#   whichNotEnd <- read.csv("whichNotEnd.csv")
#   HyperPAssign <- read.csv("HyperPAssign.csv")
#   X <- as.matrix(read.csv("X.csv"))
#   speciesToKeep <- read.csv("speciesToKeep.csv")
#   predictorDf <- read.csv("predictorDf.csv")
#
#
#   load(paste0(workPath, "/FitModel_NoEndemics.Rdata"))
#
#   #Calls coeff for predictors in matrix X and intersect
#   summ1 <-summary(FitModel, pars = c("a",paste("b[",1:dim(X)[2],"]", sep="")))$summary
#   #Calls r squareds
#   summ2 <- summary(FitModel, pars = c("r_sq","r_sq_justX","r_sq_notGeo" ))$c_summary[,,2] #[,"50%"]
#
#
#
#   a.rhat <- summ1[1,10]
#   a.mean <- summ1[1,1]
#   a.2pt5<- summ1[1,4]
#   a.97pt5 <- summ1[1,8]
#   a.sig <- (a.2pt5 * a.97pt5 > 0)
#
#   b1.rhat <- summ1[2,10]
#   b1.mean <- summ1[2,1]
#   b1.2pt5 <- summ1[2,4]
#   b1.97pt5 <- summ1[2,8]
#   b1.sig <- (b1.2pt5 * b1.97pt5 > 0)
#
#   b2.rhat <- NA
#   b2.mean <- NA
#   b2.2pt5 <-NA
#   b2.97pt5 <- NA
#   b2.sig <- NA
#
#   if (!is.na(b2)) {
#     b2.rhat <- summ1[3,10]
#     b2.mean <- summ1[3,1]
#     b2.2pt5 <- summ1[3,4]
#     b2.97pt5 <- summ1[3,8]
#     b2.sig <- (b2.2pt5 * b2.97pt5 > 0)
#   }
#
#
#   rsquared.mean <- summ2[1,1]
#   rsquared.2pt5 <- summ2[1,3]
#   rsquared.97pt5 <- summ2[1,7]
#
#
#   rsquared.justX.mean <- summ2[2,1]
#   rsquared.justX.2pt5 <- summ2[2,3]
#   rsquared.justX.97pt5 <- summ2[2,7]
#
#
#   rsquared.notGeo.mean <- summ2[3,1]
#   rsquared.notGeo.2pt5 <- summ2[3,3]
#   rsquared.notGeo.97pt5 <- summ2[3,7]
#
#   thisrow <- c(threat, b1, b2, states, scrambled, a.rhat, a.mean, a.2pt5, a.97pt5, b1.rhat, b1.mean, b1.2pt5, b1.97pt5, b1.sig, b2.rhat, b2.mean, b2.2pt5, b2.97pt5, b2.sig, rsquared.mean, rsquared.2pt5, rsquared.97pt5, rsquared.justX.mean, rsquared.justX.2pt5, rsquared.justX.97pt5, rsquared.notGeo.mean, rsquared.notGeo.2pt5, rsquared.notGeo.97pt5)
#   if (length(thisrow) < 28) return()
#   thisrow <- matrix(data = thisrow, nrow = 1, ncol = 28, dimnames = list('NA', resultnames))
#   results <- rbind(results, thisrow)
#   write.csv(results, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/results.csv')
#
#
# }
#
#
# allRoads_sriProxRoadsStd_stTRUE
# runALot('allRoads', c('sriProxRoadsStd'), TRUE, 4, FALSE)
# runALot('allRoads', c('sriProxRoadsStd'), TRUE, 4, TRUE)
#
# carDirs <- list.dirs(path = 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses', full.names = TRUE, recursive = TRUE)
# carDirs <- carDirs[ which(carDirs != 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled')]
# carDirs <- carDirs[ which(carDirs != 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/results.csv')]
# carDirs <- carDirs[ which(carDirs != 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses')]
# lapply(carDirs, readCARModel)
#
# # 6 20 17
# newCARS <- c('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c1p0x1x2x3_allDev_logCensusPopStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c1p0x1x2x3_allDev_venterPop1990sStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c1p0x1x2x3_allDev_venterPop2010sStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c2pt3_livestock_venterPastureStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/c2pt3_livestock_venterPastureStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/srambled/c2pt1_crops_cropSuitStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled/c11p3x4_tempPrecip_temAllYrsStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled/c2pt3_livestock_venterPastureStd_stTRUE', 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CountyEndangeredSpecies-master-June/CAR Analyses/scrambled/c2pt3_livestock_nlcdPastureStd_stTRUE')
# lapply(newCARS, readCARModel)
#
# #### 8gb####
# runCARModel('c2pt3_livestock', c('nlcdPastureStd'), TRUE, 4)
#
# runScrambledCARModel('c11p3x4_tempPrecip', c('kdeMeanStabStd'), TRUE, 4)
#
# runCARModel('anyThreat', c('iucn1aStd'), TRUE, 4)
#
# runCARModel('anyThreat', c('gap12X1aStd'), TRUE, 4)
#
# runCARModel('anyThreat', c('gap12Std'), TRUE, 4)
#
# runCARModel('allRoads', c('TIGERroadDensityStd'), TRUE, 4)
#
# runCARModel('allRoads', c('venterRoadStd'), TRUE, 4)
#
# runCARModel('allRoads', c('sriProxRoadsStd'), TRUE, 4)
#
# runScrambledCARModel('anyThreat', c('iucn1aStd'), TRUE, 4)
#
# runScrambledCARModel('anyThreat', c('gap12X1aStd'), TRUE, 4)
#
# runScrambledCARModel('anyThreat', c('gap12Std'), TRUE, 4)
#
# # 6 19 2017
#
# runScrambledCARModel('c2pt1_crops', c('cropSuitStd'), TRUE, 4)
#
#
#
# #### 16gb ####
#
# runScrambledCARModel('c2pt3_livestock', c('nlcdPastureStd'), TRUE, 4)
#
# runScrambledCARModel('c1p0x1x2x3_allDev', c('nlcdDevStd'), TRUE, 4)
#
# runScrambledCARModel('allRoads', c('TIGERroadDensityStd'), TRUE, 4)
#
# runScrambledCARModel('allRoads', c('venterRoadStd'), TRUE, 4)
#
# runScrambledCARModel('allRoads', c('sriProxRoadsStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('venterHFIStd'), TRUE, 2)
#
# runCARModel('c1p0x1x2x3_allDev', c('percentImpervStd'), TRUE, 2)
#
# runCARModel('c2p0x1x2x3_allAg', c('cropSuitStd'), TRUE, 2)
#
# runCARModel('c2pt1_crops', c('cropSuitStd'), TRUE, 2)
#
# runCARModel('c2pt3_livestock', c('venterPastureStd'), TRUE, 2)
#
# runCARModel('c1p0x1x2x3_allDev', c('logCensusPopStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('venterPop1990sStd'), TRUE, 4)
#
#
# #### 32 gb 2 ####
# runCARModel('c11p3x4_tempPrecip', c('temLast30Std'), TRUE, 4)
#
# runCARModel('c11p3x4_tempPrecip', c('temAllYrsStd'), TRUE, 4)
#
# # runCARModel('c1p0x1x2x3_allDev', c('logCensusPopStd', TRUE, 4)
#
# #runCARModel('c1p0x1x2x3_allDev', c('venterPop1990sStd', TRUE, 4)
#
# #runCARModel('c1p0x1x2x3_allDev', c('venterPop2010sStd', TRUE, 4)
#
# #runScrambledCARModel('c11p3x4_tempPrecip', c('tempLast30Std', TRUE, 4)
#
# runScrambledCARModel('c11p3x4_tempPrecip', c('temAllYrsStd'), TRUE, 4)
#
# #runScrambledCARModel('c1p0x1x2x3_allDev', c('logCensusPopStd', TRUE, 4)
#
# #runScrambledCARModel('c1p0x1x2x3_allDev', c('venterPop1990sStd', TRUE, 4)
#
# # 6 19 2017
# runCARModel('c1p0x1x2x3_allDev', c('venterPop2010sStd'), TRUE, 4)
#
# runScrambledCARModel('c11p3x4_tempPrecip', c('tempLast30Std'), TRUE, 4)
#
#
# #### 32 gb 3 ####
# runScrambledCARModel('c1p0x1x2x3_allDev', c('venterPop2010sStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('venterHFIStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('venterBuiltStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('venterNtLtStd'), TRUE, 4)
#
# runCARModel('c1p0x1x2x3_allDev', c('nlcdDevStd'), TRUE, 4)
#
# runScrambledCARModel('c1p0x1x2x3_allDev', c('venterHFIStd'), TRUE, 4)
#
# runScrambledCARModel('c1p0x1x2x3_allDev', c('venterBuiltStd'), TRUE, 4)
#
# runScrambledCARModel('c1p0x1x2x3_allDev', c('venterNtLtStd'), TRUE, 4)
#
# runScrambledCARModel('c1p0x1x2x3_allDev', c('nlcdDevStd'), TRUE, 4)
#
# # 6 19 2017
#
# runCARModel('c2pt3_livestock', c('venterPastureStd'), TRUE, 4)
#
# runScrambledCARModel('c2pt3_livestock', c('venterPastureStd'), TRUE, 4)
#
# #### 32 gb 4 ####
#
# runCARModel('c2p0x1x2x3_allAg', c('hanppStd'), TRUE, 4)
#
# runCARModel('c2p0x1x2x3_allAg', c('cropSuitStd'), TRUE, 4)
#
# runCARModel('c2pt1_crops', c('venterCropStd'), TRUE, 4)
#
# runCARModel('c2pt1_crops', c('nlcdCropsStd'), TRUE, 4)
#
# runCARModel('c2pt1_crops', c('cropSuitStd'), TRUE, 4)
#
# runScrambledCARModel('c2p0x1x2x3_allAg', c('hanppStd'), TRUE, 4)
#
# runScrambledCARModel('c2p0x1x2x3_allAg', c('cropSuitStd'), TRUE, 4)
#
# runScrambledCARModel('c2pt1_crops', c('venterCropStd'), TRUE, 4)
#
# runScrambledCARModel('c2pt1_crops', c('nlcdCropsStd'), TRUE, 4)
#
# # 6 19 2017
#
#
# runCARModel('c2pt3_livestock', c('nlcdPastureStd'), TRUE, 4)
#
# runScrambledCARModel('c2pt3_livestock', c('nlcdPastureStd'), TRUE, 4)
#
#



# Here is code for endemics...skip
# # tally # counties each species occurs in ####
# findSpeciesOccurrences <- function(spID){
#   spOccurrences <- grep(spID, natureServeCountyOccurrence$ELEMENT_GLOBAL_ID)
#   return(as.numeric(length(spOccurrences)))
# }
#
# speciesThreats$numberOfCounties <- sapply(speciesThreats$ELEMENT_GLOBAL_ID, findSpeciesOccurrences)
#
# # separate df for endemics ####
# speciesThreatsEndemics <- speciesThreats[ which(speciesThreats$numberOfCounties == 1), ]


# # tally number of species endemic to each county ####
# findEndemicOccurrences <- function(FIPS_CODE){
#   if (FIPS_CODE %in% natureServeCountyOccurrence$FIPS_CODE_LONG) {
#     localSpecies <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$FIPS_CODE_LONG == FIPS_CODE), ]
#     localSpecies$isEndemic <- FALSE
#     nLocal <- as.numeric(length(localSpecies$ELEMENT_GLOBAL_ID))
#     for (i in (1:nLocal)) {
#       if(localSpecies$ELEMENT_GLOBAL_ID[i] %in% speciesThreatsEndemics$ELEMENT_GLOBAL_ID) localSpecies$isEndemic[i] <- TRUE
#     }
#     nEndemic <- grep(TRUE, localSpecies$isEndemic)
#
#     return(as.numeric(length(nEndemic)))
#   }
#   else return(0)
# }
#
# natureServeCounties$endemicsInCounty <- sapply(natureServeCounties$FIPS_CODE_LONG, findEndemicOccurrences)
#
# # tally number of endemics affected by each threat in each county ####
# countEndemicsThreatened <- function(CTY_FIPS, threat){
#   natureServeCountyIndex <- as.numeric(match(CTY_FIPS, natureServeCounties$FIPS_CODE_LONG))
#   if (natureServeCounties$endemicsInCounty[natureServeCountyIndex] == 0) return(0)
#
#   speciesInCounty <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$FIPS_CODE_LONG == CTY_FIPS), ]
#   speciesInCounty <- speciesInCounty$ELEMENT_GLOBAL_ID
#
#   endemicsInCountyIndices <- match(speciesInCounty, speciesThreatsEndemics$ELEMENT_GLOBAL_ID)
#   endemicsInCountyIndices <- endemicsInCountyIndices[!is.na(endemicsInCountyIndices)]
#
#   endemicsInCounty <- speciesThreatsEndemics[endemicsInCountyIndices, ]
#
#   speciesThreatsColIndex <- as.numeric(match(threat, colnames(endemicsInCounty)))
#
#   endemicsThreatenedInCounty <- endemicsInCounty[ which(endemicsInCounty[speciesThreatsColIndex] == 1), ]
#
#   return(as.numeric(length(endemicsThreatenedInCounty$ELEMENT_GLOBAL_ID)))
# }
#
# # tally number of endemics threatened by a threat ###
# countEndemicsThreatenedByThreat <- function(threatName) {
#
#   threat <- threatName
#
#   threatCol <- apply(as.matrix(natureServeCounties$FIPS_CODE_LONG), 1, countEndemicsThreatened, threat = threat)
#
#   threatCol <- as.matrix(threatCol)
#
#   colnames(threatCol) <- paste("endemicsThreatenedBy_", threat, sep = "")
#
#   newNatureServeCounties <- cbind(natureServeCounties, threatCol)
#
#   return(newNatureServeCounties)
# }

# anyThreat <- countEndemicsThreatenedByThreat("anyThreat")
#
# natureServeCounties <- countEndemicsThreatenedByThreat("c11p3x4_tempPrecip")
#
# natureServeCounties <- countEndemicsThreatenedByThreat("c4p1p0_roads")
#
# natureServeCounties <- countSpeciesThreatenedByThreat("anyThreat")




say('+++++++++++++++++++++++++ DONE +++++++++++++++++++++++++', pre=2)
say('+++++++++++++++++++++++++ DONE +++++++++++++++++++++++++')
say('+++++++++++++++++++++++++ DONE +++++++++++++++++++++++++')
say('+++++++++++++++++++++++++ DONE +++++++++++++++++++++++++')
say('+++++++++++++++++++++++++ DONE +++++++++++++++++++++++++')
