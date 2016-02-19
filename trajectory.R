rm(list = ls())
#setwd("./Google_Drive/Research/Manuscript/2016_BNLearn/code")

##########################################################
#          function definition
##########################################################
comparePlot <- function(rawSDP, rawRD, type = "correlation", filename = "temp") {
  library(ggplot2)
  dat.SDP <- data.frame(iteration = 1:ncol(rawSDP),
                        mean = colMeans(rawSDP),
                        sd = apply(rawSDP, 2, sd))
  
  dat.RD <- data.frame(iteration = 1:ncol(rawRD),
                        mean = colMeans(rawRD),
                        sd = apply(rawRD, 2, sd))
  
  ggplot() +
        geom_line(data = dat.SDP, aes(x=iteration, y=mean, col = "SDP")) +
        geom_ribbon(data = dat.SDP, aes(x=iteration, ymin=mean-1.96*sd, ymax=mean+1.96*sd), fill = "brown1", alpha = 0.3) +
        geom_line(data = dat.RD, aes(x=iteration, y=mean, col = "RD")) +
        geom_ribbon(data = dat.RD, aes(x=iteration, ymin=mean-1.96*sd, ymax=mean+1.96*sd), fill = "cyan3", alpha = 0.3) +
        ylab(label = type) +
        theme_bw() + theme(axis.text=element_text(size=20), axis.title=element_text(size=20,face = "bold")) +
        theme(legend.position = c(0.1,0.9), legend.key = element_rect(colour = "white"), legend.text=element_text(size=20)) +
        scale_colour_manual("", 
                            breaks = c("RD", "SDP"),
                            values = c("blue", "red")) 
  ggsave(filename = paste(filename, "_", type, ".png", sep = "")) 
}

# settings
bn <- "insurance"
numExpert <- 4
sigma2 <- 1
filename <- paste('./result/', bn, '_', numExpert, '_', sigma2, sep = '')

##########################################################
#          correlation comparison
##########################################################
corr.SDP <- read.table(file = paste(filename, '_corr_SDP.csv', sep = ''), sep = ',')
corr.RD <- read.table(file = paste(filename, '_corr_RD.csv', sep = ''), sep = ',')

comparePlot(corr.SDP, corr.RD, type = "correlation", filename)

##########################################################
#          variance comparison
##########################################################
var.SDP <- read.table(file = paste(filename, '_var_SDP.csv', sep = ''), sep = ',')
var.RD <- read.table(file = paste(filename, '_var_RD.csv', sep = ''), sep = ',')

comparePlot(var.SDP, var.RD, type = "mean of variance", filename)
