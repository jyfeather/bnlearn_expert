rm(list = ls())

library(gdata)
originalData <- read.xls("/Volumes/mini/dataset/KPI/HPM_limei.xlsx", sheet = 1, header = T)

originalData[originalData=="#NULL!"] <- NA

constructs <- c(14,3,
                34,6,
                40,6,
                51,6,
                157,4,
                161,5,
                175,4,
                186,4,
                190,5,
                196,5,
                201,1,
                202,1,
                203,1,
                204,1,
                205,1,
                206,1,
                207,1,
                208,1,
                209,1,
                210,6,
                216,1,
                311,4,
                315,5,
                320,7,
                327,4,
                1429,5,
                1438,2,
                1440,2,
                1442,2,
                1444,2,
                1508,2)
constructs <- matrix(constructs, ncol = 2, byrow = T)

construct.PCA <- matrix(0, nrow = nrow(originalData), ncol = nrow(constructs))
for (i in 1:nrow(constructs)) {
  construct.start <- constructs[i,1]
  construct.end <- construct.start + constructs[i,2] - 1
  construct.value <- originalData[,construct.start:construct.end]
  
  if (construct.start == construct.end) {
    construct.value <- as.numeric(construct.value)
    construct.value[is.na(construct.value)] <- mean(construct.value, na.rm = T)
    construct.PCA[,i] <- as.numeric(construct.value)
  } else {
    construct.value <- data.matrix(construct.value)
    for (j in 1:constructs[i,2]) {
      valueset <- construct.value[,j] 
      valueset[is.na(valueset)] <- mean(valueset, na.rm = T)
      construct.value[,j] <- valueset
    }
    
    # run PCA for construct with more than 1 variable
    prc <- prcomp(data.matrix(construct.value), retx = T, center = T, scale. = T)
    res.pca <- scale(construct.value) %*% varimax(prc$rotation)$loadings
    construct.PCA[,i] <- as.numeric(res.pca[,1])
  }
}

# write
write.table(data.matrix(construct.PCA), file = "./Google_Drive/Research/Manuscript/2016_BNLearn/data/KPI.csv", 
            row.names = F, col.names = F, sep = ",")
