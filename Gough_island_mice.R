setwd("G:/Dropbox/Academic/Projects & Papers/Gough Island mice/Data")
setwd("C:/Users/Anthony Caravaggi/Dropbox/Academic/Projects & Papers/Gough Island mice/Data - revised")

library(plyr)
library(ggplot2)
library(cowplot)

bDAT = read.csv("giDAT.csv", header = TRUE, stringsAsFactors=FALSE)
bDAT["gough.se"] <- bDAT$gough.sd/sqrt(bDAT$gough.count)
bDAT["typical.se"] <- bDAT$typical.sd/sqrt(bDAT$typical.count)
bDAT["bs.diff"] <- bDAT$typical.bs-bDAT$gough.bs
bDAT["se.diff"] <- abs(bDAT$typical.se-bDAT$gough.se)

# Function to calculate minimum and maximum variance based on n randomly
# generated integers from a random distribution (nested tnorm), with truncation
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
gDAT[9,c(6:8)] <- c(150, 110, 190) # Southern giant petrel; rounded to 10, a-priori
gDAT[10,c(6:8)] <- c(550, 400, 700) # Tristan albatross; rounded to 50, a-priori
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
tDAT[9,c(6:8)] <- c(180, 150, 200) # Southern giant petrel; rounded to 10, a-priori
tDAT[10,c(6:8)] <- c(1250, 1150, 1400) # Tristan albatross; rounded to 50, a-priori
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
pDAT[9,c(6:8)] <- c(30, 0, 50) # Southern giant petrel; rounded to 10, a-priori
pDAT[10,c(6:8)] <- c(700, 650, 750) # Tristan albatross; rounded to 50, a-priori
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
sDAT["size"] <- c(520, 1800, 160, 830, 1000, 160, 250, 2100, 3800, 6800) 
sDAT$species <- factor(sDAT$species, levels = sDAT$species[order(sDAT$size)])

# Split data by each category

season <- split(sDAT, sDAT$season)
summer <- season$S
winter <- season$W
nest <- split(sDAT, sDAT$nest)
above <- nest$A
below <- nest$B

summer <- summer[order(summer$size),]
winter <- winter[order(winter$size),]
above <- above[order(above$size),]
below <- below[order(below$size),]

# t-tests between

t.test(sDAT$gough.bs~sDAT$nest,data=sDAT) 
t.test(sDAT$gough.bs~sDAT$season,data=sDAT) 

# remove outliers in great shearwater and Tristan albatross

s2DAT <- sDAT[!(sDAT$abbr=="GRSH"),]
s2DAT <- s2DAT[!(s2DAT$abbr=="TRAL"),]

t.test(s2DAT$gough.bs~s2DAT$nest,data=s2DAT) 

# Remove both prions

s3DAT <- sDAT[!(sDAT$abbr=="MAPR"),]
s3DAT <- s3DAT[!(s3DAT$abbr=="BBPR"),]

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

nest.D <- cohens_d(above$gough.bs, below$gough.bs)
season.D <- cohens_d(summer$gough.bs, winter$gough.bs)

# Season minus MacGillivray's prion
#
#season2 <- split(s2DAT, s2DAT$season)
#summer2 <- season2$S
#winter2 <- season2$W
#
# Season minus MacGillivray's and broad-billed prions
#
#season3 <- split(s3DAT, s3DAT$season)
#summer3 <- season3$S
#winter3 <- season3$W
#
#season.MAPR.D <- cohens_d(summer2$gough.bs, winter2$gough.bs)
#season.BBPR.D <- cohens_d(summer3$gough.bs, winter3$gough.bs)

# Plot  and pair plots for each category (above/below ground; summer/winter)
# Requires theme_ac1 https://github.com/arcaravaggi/browncoat/blob/master/theme_ac1.R
#
# Function for customised x axis labels (with breaks)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}


p1 <- ggplot(above, aes(x = above$species, y = above$gough.bs))+
  geom_bar(stat = "identity", width=0.85, fill = "grey", colour="black")+
  geom_errorbar(aes(ymin = above$ymin, ymax = above$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  scale_x_discrete(breaks=unique(above$species), 
                   labels=addline_format(above$species)) +
  theme_ac1() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="serif"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14,face="bold"),
        panel.spacing.y=unit(3, "lines"))+
  labs(x = "Average breeding success",y = "Species")


p2 <- ggplot(below, aes(x = below$species, y = below$gough.bs))+
  geom_bar(stat = "identity", width=0.85, fill = "grey", colour="black")+
  geom_errorbar(aes(ymin = below$ymin, ymax = below$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  scale_x_discrete(breaks=unique(below$species), 
                   labels=addline_format(below$species)) +
  theme_ac1() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="serif"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        panel.spacing.y=unit(3, "lines"))+
  labs(x = "Average breeding success",y = "Species")

p12 <- ggdraw()+
  draw_plot(p1, x = 0, y = 0.5, width = 0.65, height = 0.5) +
  draw_plot(p2, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("a)", "b)"), size = 14, x = c(0,0), y = c(1,0.5))

ggsave("plot_Nesters.png", p12, height = 10, width = 10, dpi = 600)

p3 <- ggplot(winter, aes(x = winter$species, y = winter$gough.bs))+
  geom_bar(stat = "identity", width=0.85, fill = "grey", colour="black")+
  geom_errorbar(aes(ymin = winter$ymin, ymax = winter$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  scale_x_discrete(breaks=unique(winter$species), 
                   labels=addline_format(winter$species)) +
  theme_ac1() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="serif"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        panel.spacing.y=unit(3, "lines"))+
  labs(x = "",y = "Species")

p4 <- ggplot(summer, aes(x = summer$species, y = summer$gough.bs))+
  geom_bar(stat = "identity", width=0.85, fill = "grey", colour="black")+
  geom_errorbar(aes(ymin = summer$ymin, ymax = summer$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  scale_x_discrete(breaks=unique(summer$species), 
                   labels=addline_format(summer$species)) +
  theme_ac1() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="serif"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size=14,face="bold"),
        axis.title.y = element_text(size=14,face="bold"),
        panel.spacing.y=unit(3, "lines"))+
  labs(x = "Average breeding success",y = "Species") 

p34 <- ggdraw()+
draw_plot(p3, x = 0, y = 0.5, width = 0.45, height = 0.5) +
  draw_plot(p4, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("a)", "b)"), size = 14, x = c(0,0), y = c(1,0.5))

ggsave("plot_Seasons.png", p34, height = 10, width = 10, dpi = 600)

# Horizontal bars

h1 <- ggplot(sDAT, aes(x = sDAT$species, y = sDAT$gough.bs, fill = sDAT$nest))+
  geom_bar(stat = "identity", width=0.85, colour="black")+
  coord_flip()+
  geom_errorbar(aes(ymin = sDAT$ymin, ymax = sDAT$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.7))+
  scale_x_discrete(breaks=unique(sDAT$species), 
                   labels=sDAT$species) +
  scale_fill_manual(values=c("#8c96c6", "#edf8fb"),
                    name="Nesting\nlocation",
                    breaks=c("A", "B"),
                    labels=c("Above-ground", "Below-ground"))+
  theme_ac1() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text=element_text(family="serif"),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size=16,face="bold"),
        axis.title.y = element_text(size=16,face="bold"),
        panel.spacing.y=unit(3, "lines"),
        legend.justification=c(1,0), 
        legend.position=c(1,0),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size =16),
        legend.background = element_rect(colour = 'transparent', fill = 'transparent'))+
  labs(x = "Species",y = "Mean breeding success") +
  geom_point(x = 4, y = 0.4, color = "black", shape=23, fill="black", size=5) +
  geom_point(x = 6, y = 0.475, color = "black", shape=23, fill="black", size=5) +
  geom_point(x = 10, y = 0.37, color = "black", shape=23, fill="black", size=5)


ggsave("plot_NS.png", h1, height = 7, width = 7, dpi = 600)



