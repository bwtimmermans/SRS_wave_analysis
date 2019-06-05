# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI_hist_plot.R")

# Function to illustrate confidence intervals on quatile sampling.
# This function samples from an asymmetric t-distribution and estimates a (95%) confidence interval about a chosen quantile.
# Requires quantile.CI (available).

# Inputs: n = sample size
#         q = quantile (0..1)

   source("./quantile_CI.R")
   func_plot_qq_CI <- function(n,q) {
      N <- c(100,1000,100000)
      #X11()
      par(mfrow=c(3,1),oma=c(2,7,2,2),mar=c(9,12,9,6),mgp=c(8,3,0))
      for (i in 1:3) {
         n <- N[i]
         AA <- rt(n=n,df=10,ncp=3)
         hist(AA,xlim=c(0,12),breaks=50,main=paste(n," samples from t(10,3) distribution."),xlab="",cex.main=5,cex.axis=5,cex.lab=5)
         abline(v=quantile(AA,probs=q),lwd=5,col="red")
         abline(v=sort(AA)[quantile.CI(n=n,q=q)$Interval],lwd=5,lty=3,col="red") }
      }

   png(filename = "./t-dist_samples.png", width = 2400, height = 3000)
   func_plot_qq_CI(1000,0.99)
   dev.off()

