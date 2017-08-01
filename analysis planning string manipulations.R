setwd("~/GitHub/Threat-mapping")

plan <- read.csv('analysis-planning.csv')
plan$threat <- NA
plan$predictor <- NA
plan$states <- NA
plan$CAR <- NA
plan$scramble <- NA
plan$justEndemics <- NA
plan$islands <- NA
plan$status <- NA

for (i in 1:nrow(plan)) {
line <- as.character(plan$Line[i])
  lineList <- unlist(strsplit(line, split = ","))
  plan$threat[i] <- unlist(strsplit(lineList[1], split = "="))[2]
plan$predictor[i] <- unlist(strsplit(lineList[2], split = "="))[2]
plan$states[i] <- unlist(strsplit(lineList[3], split = "="))[2]
plan$CAR[i] <- unlist(strsplit(lineList[5], split = "="))[2]
plan$scramble[i] <- unlist(strsplit(lineList[6], split = "="))[2]
plan$justEndemics[i] <- unlist(strsplit(lineList[7], split = "="))[2]
plan$islands[i] <- unlist(strsplit(lineList[5], split = "="))[2]
plan$islands[i] <- unlist(strsplit(islands, split = ")"))[1]
}

write.csv(plan, 'sofar.csv')
