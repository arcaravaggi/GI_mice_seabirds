# Function to calculate observed, minumum and maximum number of chicks
# based on boostrapped breeding success data (minmax)
# 
# dat = dataframe
# pe = argument ("obs", "min.c", "max.c")

popfunc <- function(dat, pe){
  vec <- length(dat[,1])
  pop <- dat[,2]
  ci <- dat[,3]
  min.e <- dat[,5]
  max.e <- dat[,6]
  for (i in seq_along(dat[,1])) {
    if (pe == "min.c"){
      vec[i] <- min.e[i] * pop[i]
    }
    else if (pe == "max.c"){
      vec[i] <- max.e[i] * pop[i]
    }
    else if (pe == "obs.c"){
      vec[i] <- ci[i] * pop[i]
    }
  }
  return(vec)
}

popfunc (dat = rDAT, pe = "obs.c")
popfunc (dat = rDAT, pe = "min.c")
popfunc (dat = rDAT, pe = "max.c")