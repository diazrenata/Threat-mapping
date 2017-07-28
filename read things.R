
setwd('H:/Global Change Program/Research/Multi-Threat Assessment/Analysis - Threat Mapping/CAR beginning 07132017')
source('Read stan function.R')

# 8gb
# readStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE, justEndemics = FALSE, islands = FALSE)

# 16gb
# readStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE, justEndemics = FALSE, islands = FALSE)
# readStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = TRUE, ncores = 4, CAR = FALSE, scramble = FALSE, justEndemics = FALSE, islands = TRUE)
# 32 1
#readStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE, justEndemics = FALSE, islands = FALSE)
readStanModel(threat = 'c2pt3_livestock', predictors = c('nlcdPastureStd'), states = FALSE, ncores = 4, CAR = FALSE, scramble = FALSE, justEndemics = FALSE, islands = TRUE)
readStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = TRUE, scramble = FALSE, justEndemics = FALSE, islands = FALSE)
# 32 3
readStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = F, scramble = FALSE, justEndemics = FALSE, islands = FALSE)
# readStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = TRUE, ncores = 4, CAR = F, scramble = FALSE, justEndemics = FALSE, islands = T)
# 32 4 
readStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = F, ncores = 4, CAR = F, scramble = FALSE, justEndemics = FALSE, islands = FALSE)
# readStanModel(threat = 'anyThreat', predictors = c('venterHFIStd'), states = F, ncores = 4, CAR = F, scramble = FALSE, justEndemics = FALSE, islands = T)
