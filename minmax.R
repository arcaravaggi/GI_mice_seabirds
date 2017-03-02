# Function to calculate minimum and maximum variance based on n bootstrapping 
# iterations (tnorm), with truncation
# 
# dat = dataframe
# n = number of bootstrapping iterations
# m = mean breeding success rate (proportional)
# s = SD/SE/CI (see cint function)
# l = lower bounds of truncation (default = -100)
# u = upper bounds of truncation (default = 100)
# r = round to x digits
# fn = variance output ("min" or "max")

minmax <- function(dat, n = 1000, m, s, l = -100, u = 100, r = 5, fn){
  vec <- length(dat[,1])
  mt <- dat[,3]
  sd <- dat[,4]
  for (i in seq_along(dat[,1])) {
    m <- mt[i]
    s <- sd[i]
    if (fn == "min.e"){
      vec[i] <- min(tnorm(n, m, s, l, u, r))
    }
    else if (fn == "max.e"){
      vec[i] <- max(tnorm(n, m, s, l, u, r))
    }
  }
  return(vec)
}
