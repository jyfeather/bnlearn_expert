rm(list = ls())
setwd("./Google_Drive/Research/Manuscript/2016_BNLearn/codes/")

name <- "mildew"

library(bnlearn)

load(paste("./DAGs/", name, ".rda", sep = ""))   # read bn strucutre in rda format

nodeNames <- names(bn)
adjMatrix <- amat(bn)

library(R.matlab)
writeMat(paste("~/Documents/MATLAB/DAGlearn/data/", name, ".mat", sep = ""), 
         adjMatrix = adjMatrix, nodeNames = nodeNames)
