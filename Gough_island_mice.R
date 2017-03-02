setwd("G:/Dropbox/Academic/Projects & Papers/Gough Island mice/Data")

library(plyr)

giDAT = read.csv("giDAT.csv", header = TRUE)
giDAT$abbr <- NULL

# Function for calculation of Confidence Intervals with varying critical values
# Includes option for Finite Population Correction
#
# dat = dataframe
# loc = location ("inv", invaded island(s); "non", non-invaded islands)
# z = critical value (default 1.96 [95% CI])
# n = sample size (default = 1)
# p = population size
# f = Finite Population Correction (FPC)
cint <- function(dat, loc, z = 1.96, n = 1, p, fpc){
  vec <- length(dat[,1])
  pop <- dat[,2]
  g <- dat[,3]
  t <- dat[,4]
  nt <- dat[,6]
  for(i in 1:nrow(dat)){
    n <- nt
    p <- pop
    s <- g
    f <- (p-n)/(p-1)
    if (loc == "inv"){
      if(fpc == "TRUE"){ 
        ci.fpc <- (z*(s/sqrt(n)))*sqrt(f)
        vec <- ci.fpc
      }  else {
        ci.std <- z*(s/sqrt(n))
        vec <- ci.std
      }
    }
    else if (loc == "non"){
      if(fpc == "TRUE"){   
        ci.fpc <- (z*(s/sqrt(n)))*sqrt(f)
        vec <- ci.fpc
      }  else {
        ci.std <- z*(s/sqrt(n))
        vec <- ci.std
      }
    }
  }
  print(vec)
}

# Function to calculate minimum and maximum variance based on n bootstrapping 
# iterations (nested tnorm), with truncation
# 
# dat = dataframe
# n = number of bootstrapping iterations
# m = mean breeding success rate (proportional)
# s = SD/SE/CI (see cint function)
# l = lower bounds of truncation (default = -100)
# u = upper bounds of truncation (default = 100)
# r = round to x digits
# fn = variance output ("min.e" or "max.e")
minmax <- function(dat, n = 1000, m, s, l = -100, u = 100, r = 5, fn){
  tnorm <- function(n, m, s, l, u, r) {   
    tdist <- round(rnorm(n, m, s), r)
    tdist[tdist < l] <- l
    tdist[tdist > u] <- u
    tdist
  }
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


# Add 50% confidence intervals to main dataframe

giDAT["ci.g"] <- cint(dat = giDAT, loc = "inv", z = 0.3085, fpc = "FALSE")
giDAT["ci.t"] <- cint(dat = giDAT, loc = "non", z = 0.3085, fpc = "TRUE")

# Create a data frame and calculate  estimates for Gough Island

rDAT <- data.frame(giDAT$species, giDAT$est_pop, giDAT$gough.bs, giDAT$ci.g)
rDAT <- rename(rDAT, c("giDAT.species"="species", "giDAT.est_pop" = "pop",
                       "giDAT.gough.bs"="gough.bs", "giDAT.ci.g"="ci"))

rDAT["min_err"] <- minmax(dat = rDAT, l=0, fn="min.e")
rDAT["max_err"] <- minmax(dat = rDAT, l=0, fn="max.e")
rDAT["obs.c"] <- popfunc (dat = rDAT, pe = "obs.c")
rDAT["min.c"] <- popfunc (dat = rDAT, pe = "min.c")
rDAT["max.c"] <- popfunc (dat = rDAT, pe = "max.c")


# Typical data (uninvaded islands)

tDAT <- data.frame(giDAT$species, giDAT$est_pop, giDAT$typical.bs, giDAT$ci.t)
tDAT <- rename(tDAT, c("giDAT.species"="species", "giDAT.est_pop" = "pop",
                       "giDAT.typical.bs"="typical.bs", "giDAT.ci.t"="ci"))

tDAT["min_err"] <- minmax(dat = tDAT, l=0, fn="min.e")
tDAT["max_err"] <- minmax(dat = tDAT, l=0, fn="max.e")
tDAT["obs.c"] <- popfunc (dat = tDAT, pe = "obs.c")
tDAT["min.c"] <- popfunc (dat = tDAT, pe = "min.c")
tDAT["max.c"] <- popfunc (dat = tDAT, pe = "max.c")

#Calculate difference between population estimates and relative impact of mice

giDAT["bs_diff"] <- giDAT$typical.bs-giDAT$gough.bs
giDAT["mouse.effect"] <- giDAT$est_pop*giDAT$bs_diff #estimate mouse effect (difference between estimates)

write.csv(rDAT, file = "d_gough.csv")
write.csv(tDAT, file = "d_typical.csv")
write.csv(giDAT, file = "d_population.csv")
