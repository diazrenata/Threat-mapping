
# Function to run stan model with or without CAR.
# threat = name of threat column
# predictors = vector of names of predictor columns
# ncores = integer how many cores to use
# CAR = true/false use CAR?
# scramble = true/false randomly reassign attributes to counties
# justEndemics = true/false restrict to species found in only 1 county?

runStanModel <- function(threat, predictors, states, ncores, CAR, scramble, justEndemics, islands){
  # load packages
  # load data saved from Master Analysis script
  # create inputs for stan
  # create directory for this run, save inputs
  # load and run compiled stan model
  # save the model
  # add modeled parameters to counties shapefile and save plots of spatial parameters (p_est, p_calcs, X)
  # save shapefile
  # save plots of descriptive parameters (rsq, alpha_tau, ps, likelihood, logloss)
  # save results to results.csv
  # clear workspace and quit this instance of R


  # load packages
  library(sp)
  library(rgdal)
  library(raster)
  options(java.parameters='-Xmx1g' )
  library(dismo)
  library(rgeos)
  library(geosphere)
  library(scales)
  library(RColorBrewer)
  library(doBy)
  library(phangorn)
  library(vegan)
  library(plyr)
  library(shinystan)
  library(RPushbullet)
  library(jsonlite)
  library(slam)
  library(rstan)

  # load data saved from Master Analysis script:
  # # speciesThreats with conglomerate threat columns (anyThreat, allRoads, etc) added
  # # natureServeCounties with predictors, X_zeros, and # of species threatened by each threat columns
  # # natureServeCountyOccurrence
  # # Format it correctly and filter to plants

  prefix <- function(x, len, pre='0') {

    # x		value to add leading characters (will be coerced to character class)
    # len	desired number of characters in x... will not be made shorter if x is longer than this
    # pre	value to pre-pend to x, will be repeated until nchar(x)==len

    x <- as.character(x)
    size <- nchar(x)
    if (nchar(x) < len) {
      addTo <- paste(rep(pre, each=len - size), collapse='')
      x <- paste0(addTo, x)
    }
    return(x)
  }
  speciesThreats <- read.csv("H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713/speciesThreats.csv")
  speciesThreats <- speciesThreats[,2:102]

  natureServeCountyOccurrence <- read.csv('//mbgcl02fs/usersdatavol2/CCSD/shared/Global Change Program/Research/Multi-Threat Assessment/Threatened Species Data (NatureServe)/Data/Working/00_NS_mv_CTY_G12ESA_list_201403 - County Occurrences.csv')

  # format FIPS codes correctly
  natureServeCountyOccurrence$FIPS_CODE_LONG <- apply(as.matrix(natureServeCountyOccurrence$FIPS_CD), 1, prefix, len = 5, pre = "0")

  # filter to plants
  natureServeCountyOccurrence <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$INFORMAL_TAX == "Ferns and relatives" | natureServeCountyOccurrence$INFORMAL_TAX == "Conifers and relatives" | natureServeCountyOccurrence$INFORMAL_TAX == "Hornworts" | natureServeCountyOccurrence$INFORMAL_TAX == "Liverworts" | natureServeCountyOccurrence$INFORMAL_TAX == "Mosses" | natureServeCountyOccurrence$INFORMAL_TAX == "Flowering Plants"), ]

  natureServeCounties <- shapefile('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713/natureServeCounties.shp')


  cols <- read.table('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/Data for stan/natureServeCounties0713/cols.txt')
  cols <- as.vector(cols[,1])

  colnames(natureServeCounties@data) <- cols

  natureServeCounties <- natureServeCounties[ which(natureServeCounties$STATE_NAME != 'Alaska'), ]
  natureServeCounties <- natureServeCounties[ which(natureServeCounties$STATE_NAME != 'Hawaii'), ]
  natureServeCounties <- natureServeCounties[ which(natureServeCounties$STATE_NAME != 'Puerto Rico'), ]


  # create inputs for stan

  colsToUse <- c('FIPS_CODE_LONG', predictors, 'STATE_NAME')
  inputShape <- natureServeCounties[ , colsToUse]
  colnames(inputShape@data) <- c('FIPS_CODE_LONG', 'X','STATE_NAME')
  for (i in 1:length(predictors)) {
    thisCol <- inputShape@data[ , i + 1]
    inputShape <- inputShape[which(!is.na(thisCol)), ]
  }

  states.num <- unique(inputShape$STATE_NAME)
  states.num <- as.data.frame(states.num)
  states.num$state.hyperP <- c(1:nrow(states.num))
  colnames(states.num) <- c('STATE_NAME', "STATEFP")
  inputShape@data <- join(inputShape@data, states.num, by = "STATE_NAME")



  nCounties <- nrow(inputShape)

  predictorCols <- 2:(length(predictors) + 1)


  allSpecies <- unique(speciesThreats$ELEMENT_GLOBAL_ID)

  occMatrix <- matrix(data=0, nrow = length(allSpecies), ncol = nrow(inputShape))
  for (i in 1:nrow(occMatrix)) {
    thisSpecies <- allSpecies[i]
    speciesCounties <- natureServeCountyOccurrence[ which(natureServeCountyOccurrence$ELEMENT_GLOBAL_ID == thisSpecies), 'FIPS_CODE_LONG']
    for (j in 1:length(speciesCounties)){
      countyIndex <- as.numeric(match(speciesCounties[j], inputShape$FIPS_CODE_LONG))
      occMatrix[i,countyIndex] <- 1
    }
  }

  speciesZeros <- vector(length=length(allSpecies))
  for (i in 1:length(allSpecies)) {
    speciesZeros[i] <- sum(occMatrix[i, ])
  }

  speciesToKeep <- cbind(allSpecies, speciesZeros)
  speciesToKeep <- as.data.frame(speciesToKeep)
  speciesToKeep$indices <- row.names(speciesToKeep)
  speciesToKeep <- speciesToKeep[ which(speciesToKeep$speciesZeros > 0), ]

  speciesIndices <- as.vector(as.numeric(speciesToKeep$indices))
  nSpecies <- length(speciesIndices)

  occMatrix <- occMatrix[speciesIndices, ]
  inputShape$Richness <- NA
  for (i in 1:nCounties) {
    inputShape$Richness[i] <- sum(occMatrix[,i])
  }

  SppCols <- vector(length = length(speciesIndices))

  for (i in 1:length(speciesIndices)) {
    occVect <- as.vector(occMatrix[i, ])
    inputShape@data <- cbind(inputShape@data, occVect)
    colnames(inputShape@data) <- c(colnames(inputShape@data[1:ncol(inputShape@data) - 1]), paste0('spp', speciesIndices[i]))
    SppCols[i] <- paste0('spp', speciesIndices[i])

  }

  inputShape$Richness <- rowSums(inputShape@data[,SppCols])

  threatIndex <- as.numeric(match(threat, colnames(speciesThreats)))

  EndSpec <- speciesToKeep
  EndSpec$end <- 0
  for(i in 1:nrow(speciesToKeep)){
    thisSpec <- speciesThreats[ which(speciesThreats$ELEMENT_GLOBAL_ID == speciesToKeep$allSpecies[i]), ]
    if (thisSpec[1,threatIndex]) EndSpec$end[i] <- 1
  }

  row.names(EndSpec) <- paste0('spp', EndSpec[,"indices"])

  #EndSpec <- as.matrix(EndSpec$end)

  if (states == FALSE) {
    nHyperP = 1
    HyperPAssign = vector( length = nCounties)


    for(i in 1:nCounties) {

      HyperPAssign[i] <- 1
    }
  }

  if (states == TRUE) {
    nHyperP = as.numeric(length(unique(inputShape$STATEFP)))
    HyperPAssign = inputShape$STATEFP
  }

  HyperPAssign <- as.numeric(HyperPAssign)


  setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master')

  source("shapeLopodData.R")
  library(rgeos)
  library(raster)
  library(ggplot2)
  library(slam)


  if (islands == FALSE) AdMatrixData = shapeLopodData(inputShape, fieldN = "Richness", fieldY = "Richness", keepFields = TRUE, Adjacency = TRUE)
  if (islands == TRUE) {
    AdMatrixData = shapeLopodData(inputShape, fieldN = "Richness", fieldY = "Richness", keepFields = TRUE, Adjacency = FALSE)

    AdMatrixData$geoInfo = list(sampledId=AdMatrixData$geoInfo$sampledId, notSampledId=data.frame(featureShape=vector(length=0), cellStan=vector(length = 0)) )
  }

  #Load Stan
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  #Load your data, N, y and X in this case are vectors without NAs of the same length (number of counties).
  #X is the predicting variable (any value) in a matrix format (so 1 or more variables can be used. If only one variable is being used it needs to be a matrix of 1 column)

  #Make sure there are objects nCounties and nSpecies with the number of species and counties respectively

  AllIncludedFeatures = rbind(AdMatrixData$geoInfo$sampledId,AdMatrixData$geoInfo$notSampledId)
  AllIncludedFeatures = AllIncludedFeatures[order(AllIncludedFeatures$cellStan),]

  #Create x predictors Matrix
  X = as.matrix(cbind(AdMatrixData$geoDataObject@data[AllIncludedFeatures$featureShape,"X"]))

  spOccMat = t(AdMatrixData$geoDataObject@data[AdMatrixData$geoInfo$sampledId$featureShape,SppCols])
  spOccMat = spOccMat[which(rowSums(spOccMat)>0),]

  EndSpec =  as.matrix(EndSpec[row.names(spOccMat),"end"])

  whichEnd = which(EndSpec == 1)
  whichNotEnd  = which(EndSpec == 0)

  nEnd = length(whichEnd)
  nNotEnd = length(whichNotEnd)

  #Field in shapefile with state
  nHyperP = length(unique(AdMatrixData$geoDataObject@data$STATEFP[AllIncludedFeatures$featureShape]))
  HyperPAssign = as.numeric(factor(as.character(AdMatrixData$geoDataObject@data$STATEFP[AllIncludedFeatures$featureShape])))

  HyperPKey = AdMatrixData$geoDataObject@data$STATE_NAME[AllIncludedFeatures$featureShape]
  HyperPKey = as.matrix(HyperPKey)
  HyperPKey = cbind(HyperPKey, HyperPAssign)
  HyperPKey = unique(HyperPKey)
  colnames(HyperPKey) <- c("STATE_NAME", "HyperPAssign")


  if (states == FALSE) {
    HyperPAssign = HyperPAssign / HyperPAssign
  }

  if(scramble == TRUE) {
    spOccMat = spOccMat[,sample(1:dim(spOccMat)[2], size = dim(spOccMat)[2], replace = F)]

    HyperPAssign <- HyperPAssign[sample(1:length(HyperPAssign), size = length(HyperPAssign), replace = F)]
    X <- X[sample(1:dim(X)[1], size = dim(X)[1], replace = F), 1]
    X <- as.matrix(X)

  }

  proportionEndSp = nEnd/(nEnd+nNotEnd)
  expectedRandomLogLoss = (proportionEndSp*log(proportionEndSp)+(1-proportionEndSp)*log(1-proportionEndSp))*-1

  if (CAR == TRUE) {
    stanData = list( nCounties = length(AdMatrixData$geoInfo$sampledId$featureShape),
                     sampledId = AdMatrixData$geoInfo$sampledId$cellStan,
                     nNotSampled = length(AdMatrixData$geoInfo$notSampledId$featureShape),
                     notSampledId = AdMatrixData$geoInfo$notSampledId$cellStan,
                     n = length(AdMatrixData$geoInfo$sampledId$featureShape) + length(AdMatrixData$geoInfo$notSampledId$featureShape) ,
                     W_n = dim(AdMatrixData$geoInfo$W_sparse)[1] ,
                     W_sparse = AdMatrixData$geoInfo$W_sparse ,
                     D_sparse = AdMatrixData$geoInfo$D_sparse ,
                     lambda = AdMatrixData$geoInfo$lambda_sparse  ,
                     nSpecies = nEnd+nNotEnd,
                     nEnd = nEnd,
                     nNotEnd = nNotEnd,
                     spOccMat = spOccMat,#Species by county Matrix (1 and 0)
                     endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
                     notEndSpp = whichNotEnd, #Which species are NOT endangered (numbers matching spOccMat)
                     nHyperP = nHyperP, #Number of categories (States)
                     HyperPAssign = HyperPAssign, #Vector of length counties assigning it to a state (numeric)
                     K = dim(X)[2],
                     x_pred = X,
                     proportionEndSp = proportionEndSp

    )

  }

  if (CAR == FALSE) {
    if (islands == TRUE) {
      stanData = list( nCounties = length(AdMatrixData$geoInfo$sampledId$featureShape),
                       sampledId = AdMatrixData$geoInfo$sampledId$cellStan,
                       nNotSampled = length(AdMatrixData$geoInfo$notSampledId$featureShape),
                       notSampledId = AdMatrixData$geoInfo$notSampledId$cellStan,
                       n = length(AdMatrixData$geoInfo$sampledId$featureShape) + length(AdMatrixData$geoInfo$notSampledId$featureShape) ,
                       nSpecies = nEnd+nNotEnd,
                       nEnd = nEnd,
                       nNotEnd = nNotEnd,
                       spOccMat = spOccMat,#Species by county Matrix (1 and 0)
                       endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
                       notEndSpp = whichNotEnd, #Which species are NOT endangered (numbers matching spOccMat)
                       nHyperP = nHyperP, #Number of categories (States)
                       HyperPAssign = HyperPAssign, #Vector of length counties assigning it to a state (numeric)
                       K = dim(X)[2],
                       x_pred = X,
                       proportionEndSp = proportionEndSp

      ) }
    if (islands == FALSE) {

      stanData = list( nCounties = length(AdMatrixData$geoInfo$sampledId$featureShape),
                       sampledId = AdMatrixData$geoInfo$sampledId$cellStan,
                       nNotSampled = length(vector(length = 0)),
                       notSampledId = vector(length = 0),
                       n = length(AdMatrixData$geoInfo$sampledId$featureShape) + length(AdMatrixData$geoInfo$notSampledId$featureShape) ,
                       nSpecies = nEnd+nNotEnd,
                       nEnd = nEnd,
                       nNotEnd = nNotEnd,
                       spOccMat = spOccMat,#Species by county Matrix (1 and 0)
                       endSpp = whichEnd , #Which species are endangered (numbers matching spOccMat)
                       notEndSpp = whichNotEnd, #Which species are NOT endangered (numbers matching spOccMat)
                       nHyperP = nHyperP, #Number of categories (States)
                       HyperPAssign = HyperPAssign, #Vector of length counties assigning it to a state (numeric)
                       K = dim(X)[2],
                       x_pred = X,
                       proportionEndSp = proportionEndSp
      )

    }

  }

  # create directory for this run, save inputs
  predictorString <- predictors[1]
  if (length(predictors) > 1) {
    for (i in 2:length(predictors)) {
      predictorString <- paste(predictorString, predictors[i], sep = "_")
    }
  }


  workPath <- paste0('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master/CAR Analyses/', threat, "_", predictorString, "_st", states, 'sc', scramble, 'C', CAR, "en", justEndemics, 'is', islands)
  if (dir.exists(workPath)) setwd(workPath)
  if (!dir.exists(workPath)) {
    dir.create(workPath)
    setwd(workPath)
  }

  saveRows <- c('nCounties', 'nSpecies', 'nEnd', 'nNotEnd', 'nHyperP', 'K')
  saveItems <- as.data.frame(saveRows)
  saveItems$Values <- c(nCounties, nSpecies, nEnd, nNotEnd, nHyperP, dim(X)[2])
  write.csv(saveItems, "nValues.csv", row.names=FALSE)
  rm(saveRows)
  rm(saveItems)

  write.csv(spOccMat, "spOccMat.csv", row.names = FALSE)
  write.table(whichEnd, "whichEnd.csv", row.names = FALSE)
  write.table(whichNotEnd, "whichNotEnd.csv", row.names = FALSE)
  write.table(HyperPAssign, "HyperPAssign.csv", row.names = FALSE)
  write.csv(X, "X.csv", row.names = FALSE)
  write.csv(speciesToKeep, "speciesToKeep.csv", row.names=FALSE)
  write.csv(AdMatrixData$geoDataObject@data, "geoData.csv", row.names = FALSE)
  write.csv(HyperPKey, 'HyperPKey.csv', row.names = FALSE)

  # load and run compiled stan model
  setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master')

  if (CAR == TRUE){

    StanModel = stan_model(file = "countyEndangered_NotEndemics_CAR.stan" )
  }
  if (CAR == FALSE){

    StanModel = stan_model(file = "countyEndangered_NotEndemics_notCAR.stan" )

  }
  #Runs the MCMC model
  FitModel = sampling(StanModel,
                      data = stanData,              # named list of data
                      chains = 4,                   # number of Markov chains
                      warmup = 2000,               # number of warmup iterations per chain
                      iter = 3000,                 # total number of iterations per chain (includes warm-up)
                      cores = ncores,                    # number of cores
                      refresh = 50                # show progress every 'refresh' iterations
  )

  # save the model
  save(FitModel, file= paste0(workPath, "/FitModel_NoEndemics.Rdata"))

  # add modeled parameters to counties shapefile and save plots of spatial parameters (p_est, p_calcs, X)
  resultsShape = AdMatrixData$geoDataObject


  resultsShape@data[AllIncludedFeatures$featureShape,"state_coeff"] = summary(FitModel, pars =paste("a_cat[",1:nHyperP,"]", sep=""))$c_summary[,,1][paste("a_cat[",HyperPAssign,"]", sep=""),"50%"]
  resultsShape@data[AdMatrixData$geoInfo$sampledId$featureShape,"p_est"] = summary(FitModel, pars = paste("p[",1:length(AdMatrixData$geoInfo$sampledId$featureShape),"]", sep=""))$c_summary[,,1][,"50%"]
  if (CAR == TRUE) resultsShape@data[AllIncludedFeatures$featureShape,"p_calc"] = summary(FitModel, pars = paste("sim_p[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,1][,"50%"]
  resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_justX"] = summary(FitModel, pars = paste("calc_p_justX[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,1][,"50%"]
  resultsShape@data[AllIncludedFeatures$featureShape,"p_calc_notGeo"] = summary(FitModel, pars = paste("calc_p_notGeo[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,1][,"50%"]
  if(CAR == TRUE) resultsShape@data[AllIncludedFeatures$featureShape,"geo_effect"] = summary(FitModel, pars = paste("geo_effect[",1:length(AllIncludedFeatures$featureShape),"]", sep=""))$c_summary[,,1][,"50%"]
  if (CAR == TRUE) resultsShape@data <- resultsShape@data[ , c('FIPS_CODE_LONG', 'X', 'STATE_NAME', 'STATEFP', 'Richness', 'detections', 'FeatureID', 'state_coeff', 'p_est', 'p_calc', 'p_calc_justX', 'p_calc_notGeo' ,'geo_effect')]
  if (CAR == FALSE) resultsShape@data <- resultsShape@data[ , c('FIPS_CODE_LONG', 'X', 'STATE_NAME', 'STATEFP', 'Richness', 'detections', 'FeatureID', 'state_coeff', 'p_est','p_calc_justX', 'p_calc_notGeo')]

  setwd(workPath)

  labelat = c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
  labeltext = c('0.1', '0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1.0')

  tiff(paste0(threat, "_", predictorString, "_", 'p_est.tiff'), width = 5000, height = 3000, units = "px")
  print( spplot(resultsShape, col="black", zcol="p_est", pretty = TRUE,  main= list(label = paste0(threat, '_', predictorString, " p_est"), cex = 8), colorkey = list(
    labels = list(
      at = labelat,
      labels = labeltext,
      cex = 7
    )
  )))
  dev.off()

  if (CAR == TRUE) {
    tiff(paste0(threat, "_", predictorString, "_", 'pclc.tiff'), width = 8000, height = 5000, units = "px")
    print( spplot(resultsShape, col="black", zcol="p_calc",main= list(label = paste0(threat, '_', predictorString, " p_calc"), cex = 3), pretty = TRUE, colorkey = list(
      labels = list(
        at = labelat,
        labels = labeltext,
        cex = 7
      )
    )))
    dev.off()
  }

  tiff(paste0(threat, "_", predictorString, "_", 'pclc_nG.tiff'), width = 8000, height = 5000, units = "px")
  print(spplot(resultsShape, col="black", zcol="p_calc_notGeo",main= list(label = paste0(threat, '_', predictorString, " p_calc_notGeo"), cex = 3), pretty = TRUE, colorkey = list(
    labels = list(
      at = labelat,
      labels = labeltext,
      cex = 7
    )
  )))
  dev.off()


  tiff(paste0(threat, "_", predictorString, "_", 'pc_jX.tiff'),width = 8000, height = 5000, units = "px")
  print( spplot(resultsShape, col="black", zcol="p_calc_justX", main= list(label = paste0(threat, '_', predictorString, " p_calc_justX"), cex = 3), pretty = TRUE, colorkey = list(
    labels = list(
      at = labelat,
      labels = labeltext,
      cex = 7
    )
  )))
  dev.off()

  min <- min(resultsShape$X, na.rm = TRUE)
  tenth <- (max(resultsShape$X, na.rm = TRUE) - min)/10
  labelXat = c(0:10)
  labelXat <- labelXat * tenth
  labelXat <- labelXat + min
  labelXtext = as.character(labelXat)
  labelXtext <- unlist(lapply(labelXtext, function(x) {
    y <- substr(x, 0, 4)
    return(y)
  }))


  tiff(paste0(threat, "_", predictorString, "_", 'X.tiff'),width = 8000, height = 5000, units = "px")
  print( spplot(resultsShape, col="black", zcol="X", main= list(label = paste0(threat, '_', predictorString, " X"), cex = 3),pretty = TRUE, colorkey = list(
    labels = list(
      at = labelXat,
      labels = labelXtext,
      cex = 7
    )
  )))
  dev.off()

  # save shapefile
  shapefile(resultsShape, filename = paste0(workPath, "/resultsShape.shp"))


  # save plots of descriptive parameters (rsq, alpha_tau, ps, likelihood, logloss)
  tiff(paste0(threat, "_", predictorString, "_", 'rsq.tiff'), width = 1000, height = 1000, units = "px")
  if (CAR == TRUE) print(stan_dens(FitModel,pars = c("r_sq","r_sq_justX", "r_sq_notGeo"),separate_chains = T ) +  labs(title = paste0(threat, '_', predictorString, " rsq")))
  if (CAR == FALSE) print(stan_dens(FitModel,pars = c("r_sq_justX", "r_sq_notGeo"),separate_chains = T ))
  dev.off()

  if (CAR == TRUE){

    tiff(paste0(threat, "_", predictorString, "_", 'a_t.tiff'), width = 1000, height = 1000, units = "px")
    print(stan_dens(FitModel,pars = c("alpha","tau"),separate_chains = T )+  labs(title = paste0(threat, '_', predictorString, " alpha_tau")))
    dev.off()

  }

  tiff(paste0(threat, "_", predictorString, "_", 'ps.tiff'), width = 1000, height = 1000, units = "px")
  print(stan_dens(FitModel,pars = c("p[1]","p[50]"),separate_chains = T )+  labs(title = paste0(threat, '_', predictorString, " ps")))
  dev.off()

  tiff(paste0(threat, "_", predictorString, "_", 'lklhd.tiff'), width = 1000, height = 1000, units = "px")
  print(stan_plot(FitModel, pars = c("random_lh", "obs_lh"), show_density=T, ci_level=0.95, outer_level=1)+  labs(title = paste0(threat, '_', predictorString, " likelihood")))
  dev.off()

  tiff(paste0(threat, "_", predictorString, "_", 'lglss.tiff'), width = 1000, height = 1000, units = "px")

  if (CAR == TRUE)  print(stan_plot(FitModel, pars = c("logloss_random", "logloss_obs", "logloss_calc_p", "logloss_notGeo", "logloss_justX" ),show_density=T, ci_level=0.95, outer_level=1)+ geom_vline(xintercept=expectedRandomLogLoss, linetype="dashed", color = "red") +  labs(title = paste0(threat, '_', predictorString, " logloss")))

  if (CAR == FALSE)  print(stan_plot(FitModel, pars = c("logloss_random", "logloss_obs", "logloss_notGeo", "logloss_justX" ), show_density=T, ci_level=0.95, outer_level=1)+ geom_vline(xintercept=expectedRandomLogLoss, linetype="dashed", color = "red")+  labs(title = paste0(threat, '_', predictorString, " logloss")))

  dev.off()


  # save results to results.csv
  results <- read.csv('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master/results.csv')
  resultnames <- c('threat', 'b1', 'b2', 'states', 'scrambled', 'CAR', 'justEndemics', 'islands', 'a.rhat', 'a.mean', 'a.2pt5', 'a.97pt5', 'b1.rhat', 'b1.mean', 'b1.2pt5', 'b1.95pt5', 'b1.sig', 'b2.rhat', 'b2.mean', 'b2.2pt5', 'b2.97pt5', 'b2.sig', 'rsquared.mean', 'rsquared.2pt5', 'rsquared.97.5', 'rsquared.justX.mean', 'rsquared.justX.2pt5', 'rsquared.justX.97pt5', 'rsquared.notGeo.mean', 'rsquared.notGeo.2pt5', 'rsquared.notGeo.97.5',
                   'random_lh_mean', 'random_lh_2pt5', 'obs_lh_mean', 'obs_lh_97pt5', 'logloss_random_mean', 'logloss_random_2pt5', 'logloss_obs_mean', 'logloss_obs_97pt5', 'logloss_calc_p_mean' ,'logloss_calc_p_97pt5', 'logloss_notGeo_mean', 'logloss_notGeo_97pt5', 'logloss_justX_mean', 'logloss_justX_97pt5', 'rsquared_log.mean', 'rsquared_log.2pt5', 'rsquared_log.97.5', 'rsquared.justX_log.mean',
                   'rsquared.justX_log.2pt5', 'rsquared.justX_log.97pt5', 'rsquared.notGeo_log.mean', 'rsquared.notGeo_log.2pt5', 'rsquared.notGeo_log.97.5', 'alpha_mean')
  results <- results[,resultnames]

  summ1 <- summary(FitModel, pars = c("a",paste("b[",1:dim(X)[2],"]", sep="")))$summary #[,"50%"]
  summ4 <- summary(FitModel, pars = c("random_lh", "obs_lh"),  use_cache = F)$summary

  if (CAR == TRUE)  {
    summ2<- summary(FitModel, pars = c("r_sq","r_sq_notGeo","r_sq_justX"))$summary
    summ3 <- summary(FitModel, pars = c("r_sq_log","r_sq_notGeo_log","r_sq_justX_log"))$summary
    summ5 <- summary(FitModel, pars = c("logloss_random", "logloss_obs", "logloss_calc_p", "logloss_notGeo", "logloss_justX" ),  use_cache = F)$summary
    summ6 <- summary(FitModel, pars = 'alpha')$summary

  }

  if (CAR == FALSE) {
    summ2<- summary(FitModel, pars = c("r_sq_notGeo","r_sq_justX"))$summary
    summ3 <- summary(FitModel, pars = c("r_sq_notGeo_log","r_sq_justX_log"))$summary
    summ5 <- summary(FitModel, pars = c("logloss_random", "logloss_obs", "logloss_notGeo", "logloss_justX" ),  use_cache = F)$summary

  }

  scrambled <- scramble
  CAR <- CAR
  justEndemics <- justEndemics
  b1 <- predictors[1]
  b2 <- NA
  if (length(predictors) > 1) b2 <- predictors[2]


  a.rhat <- summ1[1,10]
  a.mean <- summ1[1,1]
  a.2pt5<- summ1[1,4]
  a.97pt5 <- summ1[1,8]
  a.sig <- (a.2pt5 * a.97pt5 > 0)

  b1.rhat <- summ1[2,10]
  b1.mean <- summ1[2,1]
  b1.2pt5 <- summ1[2,4]
  b1.97pt5 <- summ1[2,8]
  b1.sig <- (b1.2pt5 * b1.97pt5 > 0)

  b2.rhat <- NA
  b2.mean <- NA
  b2.2pt5 <-NA
  b2.97pt5 <- NA
  b2.sig <- NA

  if (!is.na(b2)) {
    b2.rhat <- summ1[3,10]
    b2.mean <- summ1[3,1]
    b2.2pt5 <- summ1[3,4]
    b2.97pt5 <- summ1[3,8]
    b2.sig <- (b2.2pt5 * b2.97pt5 > 0)
  }

  if ( CAR == TRUE) {

    rsquared.mean <- summ2[1,1]
    rsquared.2pt5 <- summ2[1,4]
    rsquared.97pt5 <- summ2[1,8]


    rsquared.justX.mean <- summ2[3,1]
    rsquared.justX.2pt5 <- summ2[3,4]
    rsquared.justX.97pt5 <- summ2[3,8]


    rsquared.notGeo.mean <- summ2[2,1]
    rsquared.notGeo.2pt5 <- summ2[2,4]
    rsquared.notGeo.97pt5 <- summ2[2,8]

    rsquared_log.mean <- summ3[1,1]
    rsquared_log.2pt5 <- summ3[1,4]
    rsquared_log.97pt5 <- summ3[1,8]


    rsquared.justX_log.mean <- summ3[3,1]
    rsquared.justX_log.2pt5 <- summ3[3,4]
    rsquared.justX_log.97pt5 <- summ3[3,8]


    rsquared.notGeo_log.mean <- summ3[2,1]
    rsquared.notGeo_log.2pt5 <- summ3[2,4]
    rsquared.notGeo_log.97pt5 <- summ3[2,8]

    logloss_random_mean <- summ5[1,1]
    logloss_random_2pt5 <- summ5[1,4]
    logloss_obs_mean <- summ5[2,1]
    logloss_obs_97pt5 <- summ5[2,8]
    logloss_calc_p_mean <- summ5[3,1]
    logloss_calc_p_97pt5 <- summ5[3,8]
    logloss_notGeo_mean <- summ5[4,1]
    logloss_notGeo_97pt5 <- summ5[4,8]
    logloss_justX_mean <- summ5[5,1]
    logloss_justX_97pt5 <- summ5[5,8]
    alpha_mean <- summ6[1,1]
  }

  if ( CAR == FALSE) {

    rsquared.justX.mean <- summ2[2, 1]
    rsquared.justX.2pt5 <- summ2[2,4]
    rsquared.justX.97pt5 <- summ2[2,8]


    rsquared.notGeo.mean <- summ2[1,1]
    rsquared.notGeo.2pt5 <- summ2[1,4]
    rsquared.notGeo.97pt5 <- summ2[1,8]

    rsquared.mean <- NA
    rsquared.2pt5 <- NA
    rsquared.97pt5 <- NA

    rsquared.notGeo_log.mean <- summ3[1,1]
    rsquared.notGeo_log.2pt5 <- summ3[1,4]
    rsquared.notGeo_log.97pt5 <- summ3[1,8]
    rsquared_log.mean <- NA
    rsquared_log.2pt5 <- NA
    rsquared_log.97pt5 <- NA

    rsquared.justX_log.mean <- summ3[2,1]
    rsquared.justX_log.2pt5 <- summ3[2,4]
    rsquared.justX_log.97pt5 <- summ3[2,8]

    logloss_random_mean <- summ5[1,1]
    logloss_random_2pt5 <- summ5[1,4]
    logloss_obs_mean <- summ5[2,1]
    logloss_obs_97pt5 <- summ5[2,8]
    logloss_calc_p_mean <- NA
    logloss_calc_p_97pt5 <- NA
    logloss_notGeo_mean <- summ5[3,1]
    logloss_notGeo_97pt5 <- summ5[3,8]
    logloss_justX_mean <- summ5[4,1]
    logloss_justX_97pt5 <- summ5[4,8]

    alpha_mean <- NA

  }

  random_lh_mean <- summ4[1,1]
  random_lh_2pt5 <- summ4[1,4]
  obs_lh_mean <- summ4[2,1]
  obs_lh_97pt5 <- summ4[2,4]


  thisrow <- c(threat, b1, b2, states, scrambled, CAR, justEndemics, islands, a.rhat, a.mean, a.2pt5, a.97pt5, b1.rhat, b1.mean, b1.2pt5, b1.97pt5, b1.sig, b2.rhat, b2.mean, b2.2pt5, b2.97pt5, b2.sig, rsquared.mean, rsquared.2pt5, rsquared.97pt5, rsquared.justX.mean, rsquared.justX.2pt5, rsquared.justX.97pt5, rsquared.notGeo.mean, rsquared.notGeo.2pt5, rsquared.notGeo.97pt5,
               random_lh_mean, random_lh_2pt5, obs_lh_mean, obs_lh_97pt5, logloss_random_mean, logloss_random_2pt5, logloss_obs_mean, logloss_obs_97pt5, logloss_calc_p_mean ,logloss_calc_p_97pt5, logloss_notGeo_mean, logloss_notGeo_97pt5, logloss_justX_mean, logloss_justX_97pt5, rsquared_log.mean, rsquared_log.2pt5, rsquared_log.97pt5, rsquared.justX_log.mean,
               rsquared.justX_log.2pt5, rsquared.justX_log.97pt5, rsquared.notGeo_log.mean, rsquared.notGeo_log.2pt5, rsquared.notGeo_log.97pt5, alpha_mean)
  if (length(thisrow) < 28) return()

  thisrow <- matrix(data = thisrow, nrow = 1, ncol = length(resultnames), dimnames = list(NULL, resultnames))
  results <- rbind(results, thisrow)
  write.csv(results, 'H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017/CountyEndangeredSpecies-master/results.csv', row.names = FALSE)
  # save states effects
  if (states == T) {
    stateeffects <- stan_plot(FitModel, pars = paste("a_cat[",1:nHyperP,"]", sep=""), show_density=T, ci_level=0.95, outer_level=1)+ geom_vline(xintercept=0, linetype="dashed", color = "red")

    HyperPKey <- read.csv("HyperPKey.csv")
    HyperPKey$parname <- paste0('a_cat[', HyperPKey$HyperPAssign, ']')
    statenames <- vector(length=length(stateeffects$data$params))
    parnames <- stateeffects$data$params


    for (i in 1:length((stateeffects$data$params))) {
      statenames[i] <- as.character(HyperPKey[ which(HyperPKey$parname == (stateeffects$data$params)[i]), 1])
      statenames[i] <- gsub(" ", "_", statenames[i])
    }
    HyperPKey$statenames <- statenames


    test <- cbind(parnames, statenames)

    # need to take the spaces out of statenames

    #levels(stateeffects$data$params) <- statenames


    allpars <- as.matrix(names(FitModel))
    row.names(allpars) <- allpars
    allpars[,1] <- 1:length(allpars)


    these <- as.numeric(allpars[parnames, ])
    names(FitModel)[these] <- statenames
    stateeffects <- stan_plot(FitModel, pars = names(FitModel)[these], show_density=T, ci_level=0.95, outer_level=1)+ geom_vline(xintercept=0, linetype="dashed", color = "red")

    tiff(paste0(threat, "_", predictorString, "_", 'steff.tiff'), width = 1000, height = 1000, units = "px")
    print(
      stateeffects
    )
    dev.off()
  }

  # create a marker file
  done <- 'done'
  write.table(done, 'done.txt')
  # clear workspace and quit this instance of R
  rm(list=ls())
  gc()
  quit(save = 'no', status = 0)
}
