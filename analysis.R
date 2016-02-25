rm(list = ls())
setwd(dir = "./Google_Drive/Research/Manuscript/2016_BNLearn/codes/")

# consider all combinations
dataNameSet <- as.factor(c('earthquake', 'asia', 'child', 'insurance', 'mildew', 'alarm', 'barley', 'hailfinder'))
numExpertSet <- as.factor(c(1, 2, 4))
sigma2Set <- as.factor(c(1, 2, 4))

combinations <- expand.grid(dataNameSet, numExpertSet, sigma2Set)
colnames(combinations) <- c("data", "numExpert", "sigma2")

# store records
resultMat <- matrix(0, nrow = nrow(combinations), ncol = 2*3)

interestIter <- c(3,6,9)
numComb <- 72
i <- 1
while (i <= numComb) {
  fileSDP <- paste("./result/", combinations[i,1], "_", combinations[i,2], "_", combinations[i,3], "_var_SDP.csv", sep = "")
  fileRD <- paste("./result/", combinations[i,1], "_", combinations[i,2], "_", combinations[i,3], "_var_RD.csv", sep = "")
  
  varSDP <- as.matrix(read.table(fileSDP, sep = ","))
  varRD <- as.matrix(read.table(fileRD, sep = ","))
  
  varSDP <- varSDP[,c(1,interestIter)]
  varRD <- varRD[,c(1,interestIter)]
  
  #diffSDP <- varSDP[,2:4] - varSDP[,1]
  #diffRD <- varRD[,2:4] - varRD[,1]
  
  #diffSDP <- t(diff(t(varSDP), lag = 1))
  #diffRD <- t(diff(t(varRD), lag = 1))
  
  #diffRatio <- diffSDP / diffRD
  
  diffTwo <- varSDP[,-1] - varRD[,-1]
  diffRD <- t(diff(t(varRD), lag = 1))
  diffRatio <- diffTwo / diffRD
  
  resultMat[i, c(1,3,5)] <- colMeans(diffRatio)
  resultMat[i, c(2,4,6)] <- apply(diffRatio, 2, sd)
  
  i <- i + 1
}

resultDF <- round(resultMat, 2)
sepDF <- matrix(100, nrow = nrow(combinations), ncol = 5)
sepDF[,c(1,3,5)] <- 50

result <- cbind(resultDF[,1],sepDF[,1],resultDF[,2],sepDF[,2],resultDF[,3],sepDF[,3],
                resultDF[,4],sepDF[,4],resultDF[,5],sepDF[,5],resultDF[,6])

write.table(result, file = "./result_analysis/diffVar.txt", sep = " ", row.names = F, col.names = F)
