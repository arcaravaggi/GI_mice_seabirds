setwd("G:/Dropbox/Academic/Projects & Papers/Gough Island mice/Data")
setwd("C:/Users/Anthony Caravaggi/Dropbox/Academic/Projects & Papers/Gough Island mice/Data")

library(plyr)
library(ggplot2)

bDAT = read.csv("giDAT.csv", header = TRUE)
bDAT["gough.se"] <- bDAT$gough.sd/sqrt(bDAT$gough.count)
bDAT["typical.se"] <- bDAT$typical.sd/sqrt(bDAT$typical.count)
bDAT["bs.diff"] <- bDAT$typical.bs-bDAT$gough.bs
bDAT["se.diff"] <- abs(bDAT$typical.se-bDAT$gough.se)

# Function to calculate minimum and maximum variance based on n randomly
# generated integers (nested tnorm), with truncation
# 
# dat = dataframe
# n = number of randomly generate integers
# m = mean breeding success rate (proportional)
# s = SD/SE/CI (see cint function)
# l = lower bounds of truncation (default = -100)
# u = upper bounds of truncation (default = 100)
# r = round to x digits
# fn = variance output ("min.e" or "max.e")
minmax <- function(dat, n = 1000, m, s, l = -100, u = 100, r = 5, fn){
  vec <- length(dat[,1])
  tnorm <- function(n, m, s, l, u, r) {   
    tdist <- round(rnorm(n, m, s), r)
    tdist[tdist < l] <- l
    tdist[tdist > u] <- u
    tdist
  }
  for (i in seq_along(dat[,1])) {
    mn <- m[i]
    sd <- s[i]
    if (fn == "min.e"){
      vec[i] <- min(tnorm(n, mn, sd, l, u, r))
    }
    else if (fn == "max.e"){
      vec[i] <- max(tnorm(n, mn, sd, l, u, r))
    }
  }
  return(vec)
}

# Function to calculate observed, minumum and maximum number of chicks
# based on normally-distributed breeding success data (minmax)
# 
# dat = dataframe
# pop = estimated population size
# bs = estimated breeding success
# r = round to nearest (default = 1000)
# min.e = minimum value of normal distribution (see minmax)
# max.e = maximum value of normal distribution (see minmax)
# pe = argument ("obs", "min.c", "max.c")
popfunc <- function(dat, pop, bs, r = 1000, min.e, max.e, pe){
  vec <- length(dat[,1])
  for (i in seq_along(dat[,1])) {
    if (pe == "min.c"){
      vec[i] <- min.e[i] * pop[i]
      vec[i] <- round(vec[i]/r)*r
    }
    else if (pe == "max.c"){
      vec[i] <- max.e[i] * pop[i]
      vec[i] <- round(vec[i]/r)*r
    }
    else if (pe == "obs.c"){
      vec[i] <- bs[i] * pop[i]
      vec[i] <- round(vec[i]/r)*r
    }
  }
  return(vec)
}

# Create a data frame and calculate  estimates for Gough Island

gDAT <- data.frame(bDAT$species, bDAT$est_pop, bDAT$gough.bs)
gDAT <- rename(gDAT, c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.gough.bs"="gough.bs"))
gDAT["min_err"] <- minmax(dat = gDAT, m = gDAT$gough.bs, s = bDAT$gough.se, l=0, fn="min.e")
gDAT["max_err"] <- minmax(dat = gDAT, m = gDAT$gough.bs, s = bDAT$gough.se, u = 1, fn="max.e")
gDAT["obs.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, pe = "obs.c")
gDAT["min.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, min.e = gDAT$min_err, pe = "min.c")
gDAT["max.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, max.e = gDAT$max_err, pe = "max.c")
gDAT["mean.c"] <- ((gDAT$max.c-gDAT$obs.c)+
                         (gDAT$obs.c-gDAT$min.c))/2

# Typical data (uninvaded islands)

tDAT <- data.frame(bDAT$species, bDAT$est_pop, bDAT$typical.bs)
tDAT <- rename(tDAT, c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.typical.bs"="typical.bs"))

tDAT["min_err"] <- minmax(dat = tDAT, m = tDAT$typical.bs, s = bDAT$typical.se, l=0, fn="min.e")
tDAT["max_err"] <- minmax(dat = tDAT, m = tDAT$typical.bs, s = bDAT$typical.se, u=1, fn="max.e")
tDAT["obs.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, pe = "obs.c")
tDAT["min.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, min.e = tDAT$min_err, pe = "min.c")
tDAT["max.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, max.e = tDAT$max_err, pe = "max.c")
tDAT["mean.c"] <- ((tDAT$max.c-tDAT$obs.c)+
                     (tDAT$obs.c-tDAT$min.c))/2


# Calculate difference between population estimates and relative impact of mice

pDAT <- data.frame(bDAT$species, bDAT$est_pop, bDAT$bs.diff)
pDAT <- rename(pDAT, c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.bs.diff"="difference.bs"))
pDAT["min_err"] <- minmax(dat = pDAT, m = pDAT$difference.bs, s = bDAT$se.diff, l=0, fn="min.e")
pDAT["max_err"] <- minmax(dat = pDAT, m = pDAT$difference.bs, s = bDAT$se.diff, u=1, fn="max.e")
pDAT["obs.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, pe = "obs.c")
pDAT["min.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, min.e = pDAT$min_err, pe = "min.c")
pDAT["max.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, max.e = pDAT$max_err, pe = "max.c")
pDAT["mean.c"] <- ((pDAT$max.c-pDAT$obs.c)+
                     (pDAT$obs.c-pDAT$min.c))/2

# Write to csv

write.csv(gDAT, file = "d_gough.csv")
write.csv(tDAT, file = "d_typical.csv")
write.csv(pDAT, file = "d_population.csv")

# Perform Chi-square test

pop <- data.frame(gDAT$obs.c, tDAT$obs.c) 
chisq.test(pop)

#####################################################
 
# Now to plot breeding success for above/below-ground breeders 
# and summer/winter breeders

sDAT<-bDAT
sDAT["ymin"] <- sDAT$gough.bs-sDAT$gough.se
sDAT["ymax"] <- sDAT$gough.bs+sDAT$gough.se

# Split data by each category

season <- split(sDAT, sDAT$season)
summer <- season$S
winter <- season$W
nest <- split(sDAT, sDAT$nest)
above <- nest$A
below <- nest$B


# t-tests between

t.test(sDAT$gough.bs~sDAT$nest,data=sDAT) 
t.test(sDAT$gough.bs~sDAT$season,data=sDAT) 

# Remove both prions

s2DAT <- sDAT[!(sDAT$abbr=="MAPR"),]
s3DAT <- s2DAT[!(s2DAT$abbr=="BBPR"),]

t.test(s2DAT$gough.bs~s2DAT$season,data=s2DAT)
t.test(s3DAT$gough.bs~s3DAT$season,data=s3DAT)


# Estimate effect size, D
# Function to calculate D
# x = x data
# y = y data

cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))       
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     
  cd  <- md/csd             
}

nest.D <- cohensD(above$gough.bs, below$gough.bs)

# Season minus MacGillivray's prion

season2 <- split(s2DAT, s2DAT$season)
summer2 <- season2$S
winter2 <- season2$W

# Season minus MacGillivray's and broad-billed prions

season3 <- split(s3DAT, s3DAT$season)
summer3 <- season3$S
winter3 <- season3$W

season.MAPR.D <- cohensD(summer2$gough.bs, winter2$gough.bs)
season.BBPR.D <- cohensD(summer3$gough.bs, winter3$gough.bs)

# Plot  and pair plots for each category (above/below ground; summer/winter)
#
# dat = dataframe
# x = categorical data
# y = continuous data
# w = bar width
# min = minimum error bar value
# max = maximum error bar value
# labx = x axis label
# laby = y axis label
# text.size.x = a axis text size
# text.size.y = y axis text size
# title.size = axis label size
#
# Annotations are used for corner labels
# xa = annotation x position
# ya = annotation y position
# ann.size = annotation text size
# laba = annotation text

bplot <- function(dat, x, y, w=0.85, min, max, labx, laby, text.size.x, 
                  text.size.y, title.size, xa, ya, ann.size, laba) {
  localenv <- environment()
  ggplot(dat, aes(x, y))+
    geom_bar(stat = "identity", width=w, fill = "grey", colour="black")+
    geom_errorbar(aes(ymin = min, ymax = max),width=w/2, position=position_dodge(.9))+
    scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
    scale_x_discrete(breaks=unique(x), 
                     labels=addline_format(x)) +
    theme_bw() + 
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          text=element_text(family="serif"),
          axis.text.x = element_text(size = text.size.x, angle = 35, hjust = 1),
          axis.text.y = element_text(size = text.size.y),
          axis.title.x = element_text(size=title.size,face="bold", hjust = 1.5),
          axis.title.y = element_text(size=title.size,face="bold"),
          panel.spacing.y=unit(3, "lines"))+
    annotate("text", x = xa, y = ya, size = ann.size, label = laba)+
    labs(x=labx,y=laby)
}

# Function for customised x axis labels (with breaks)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}



p1 <- bplot(dat = above, x = above$species, y = above$gough.bs, w = 0.85, min = above$ymin, max = above$ymax,
            labx="Species", laby = "Average breeding success", text.size.x = 14, text.size.y = 14, title.size = 16, 
            xa = 0.65, ya = 0.67, ann.size = 5, laba = "a)")

p2 <- bplot(dat = below, x = below$species, y = below$gough.bs, w = 0.85, min = below$ymin, max = below$ymax,
            labx="", laby = "", text.size.x = 14, text.size.y = 0, title.size = 16,
            xa = 0.65, ya = 0.67, ann.size = 5, laba = "b)")

p3 <- bplot(dat = summer, x = summer$species, y = summer$gough.bs, w = 0.85, min = summer$ymin, max = summer$ymax,
            labx="Species", laby = "Average breeding success", text.size.x = 14, text.size.y = 14, title.size = 16,
            xa = 0.65, ya = 0.67, ann.size = 5, laba = "a)")

p4 <- bplot(dat = winter, x = winter$species, y = winter$gough.bs, w = 0.85, min = winter$ymin, max = winter$ymax,
            labx="", laby = "", text.size.x = 14, text.size.y = 0, title.size = 16,
            xa = 0.65, ya = 0.67, ann.size = 5, laba = "b)")

# Multiple plot function
#
# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Export multiplot

png("plot_Nesters.png", width = 12, height = 4, units = 'in', res = 600)
multiplot(p1, p2, cols=2)
dev.off()

png("plot_Seasons.png", width = 12, height = 4, units = 'in', res = 300)
multiplot(p3, p4, cols=2)
dev.off()


# Could also use facet plots, as below, though adding corner labels is trickier

ggplot(sDAT, aes(sDAT$abbr, sDAT$gough.bs))+
  geom_bar(stat = "identity", width=0.85, fill = "grey", colour="black")+
  geom_errorbar(aes(ymin = sDAT$ymin, ymax = sDAT$ymax),width=0.5, position=position_dodge(.9))+
  facet_grid(. ~ sDAT$nest, scales = "free", space = "free")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"),
        text=element_text(family="serif"),
        strip.text.x = element_blank(),
        axis.text.x = element_text(size = 16, angle = 30, hjust = 1),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        panel.spacing.y=unit(3, "lines"))+
  labs(x="Species",y="Average breeding success")

