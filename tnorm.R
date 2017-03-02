# Function for truncation of normal distribution
# n = number of iterations (default = 1000)
# m = mean
# s = SD/SE/CI
# l = lower bounds (default = -100)
# u = upper bounds (default = 100)
# r = round to x integers (default = 5)
tnorm <- function(n = 1000, m, s, l = -100, u = 100, r = 5) {   
  tdist <- round(rnorm(n, m, s), r)
  tdist[tdist < l] <- l
  tdist[tdist > u] <- u
  tdist
}