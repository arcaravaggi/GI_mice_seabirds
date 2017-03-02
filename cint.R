setwd("G:/Dropbox/Academic/Projects & Papers/Gough Island mice/Data")

giDAT = read.csv("giDAT.csv", header = TRUE)

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



#cint(giDAT, z = 1.96, fpc = "TRUE")  #95%
#cint(giDAT, z = 1.28, fpc = "TRUE")  #80%
#cint(giDAT, z = 1.03, fpc = "TRUE")  #75%
#cint(giDAT, z = 0.84, fpc = "TRUE")  #70%
#cint(giDAT, z = 0.675, fpc = "TRUE")  #65%

output <- cint(dat = giDAT, loc = "inv", z = 0.3085, fpc = "TRUE")
giDAT["cint"]   <- output
