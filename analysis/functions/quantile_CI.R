# source("/home/ben/research/NOC/SRS_wave_analysis/analysis/functions/quantile_CI.R")
#
# Near-symmetric distribution-free confidence interval for a quantile `q`.
# Returns indexes into the order statistics.
#
quantile.CI <- function(n, q, alpha=0.05) {
  #
  # Search over a small range of upper and lower order statistics for the 
  # closest coverage to 1-alpha (but not less than it, if possible).
  #
  u <- qbinom(1-alpha/2, n, q) + (-2:2) + 1
  l <- qbinom(alpha/2, n, q) + (-2:2)
  u[u > n] <- Inf
  l[l < 0] <- -Inf
  coverage <- outer(l, u, function(a,b) pbinom(b-1,n,q) - pbinom(a-1,n,q))
  if (max(coverage) < 1-alpha) i <- which(coverage==max(coverage)) else
    i <- which(coverage == min(coverage[coverage >= 1-alpha]))
  i <- i[1]
  #
  # Return the order statistics and the actual coverage.
  #
  u <- rep(u, each=5)[i]
  l <- rep(l, 5)[i]
  return(list(Interval=c(l,u), Coverage=coverage[i]))
}
#
# Example: test coverage via simulation.
#
#n <- 100      # Sample size
#q <- 0.90     # Percentile
#
# You only have to compute the order statistics once for any given (n,q).
#
#lu <- quantile.CI(n, q)$Interval
#
# Generate many random samples from a known distribution and compute 
# CIs from those samples.
#
   #set.seed(17)
   #n.sim <- 1e4
   #index <- function(x, i) ifelse(i==Inf, Inf, ifelse(i==-Inf, -Inf, x[i]))
   #sim <- replicate(n.sim, index(sort(rnorm(n)), lu))
#
# Compute the proportion of those intervals that cover the percentile.
#
   #F.q <- qnorm(q)
   #covers <- sim[1, ] <= F.q & F.q <= sim[2, ]
#
# Report the result.
#
   #message("Simulation mean coverage was ", signif(mean(covers), 4), 
   #     "; expected coverage is ", signif(quantile.CI(n,q)$Coverage, 4))

