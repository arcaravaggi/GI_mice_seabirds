library(ggplot2)
library(Rmisc)
library(lme4)
library(lsmeans)
library(scales)
library(broom)
library(plyr)
library(dplyr)
library(MuMIn)

# Read data for Gough species and analogues
anDat <- read.csv("gough_analogue.csv", header = TRUE)

# Filter new df for Gough species only
gDat <- anDat %>% filter(l2 == "GI")
b <- data.frame(species = unique(gDat$species),
                abb = c("ATPE", "AYNA", "BRBP", "GRPE", "GRSH", "MGPE", "SGPE", "SOAL", "SOPE", "TRAL"))
gDat <- merge(gDat, b, by = "species")

# Summarise breeding success, maintaining associated data
gsum <- summarySE(gDat, measurevar="bs", groupvars=c("species", "abb", "weight", "season", "nest"))
ansum <- summarySE(anDat, measurevar="bs", groupvars=c("species", "l2"))

# Plot theme
theme_ac1 <- function(base_family = "serif", base_size_a = 12, base_size_t = 12){
  theme_bw(base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.grid = element_blank(),   
      axis.text = element_text(size = base_size_a),
      axis.title = element_text(size=base_size_t,face="bold"),
      legend.key=element_rect(colour=NA, fill =NA),
      panel.border = element_rect(fill = NA, colour = "black", size=0),
      panel.background = element_rect(fill = "white", colour = "black"), 
      strip.background = element_rect(fill = NA)
    )
}

# Plot Gough birds - weight*reproductive success, with nesting location and season as levels

gsum["ymin"] <- gsum$bs-gsum$se
gsum["ymax"] <- gsum$bs+gsum$se
gsum$ymin[3] <- 0 # Replace values <0 with 0
gsum$ymin[6] <- 0 # Replace values <0 with 0
gsum$species <- factor(gsum$species, levels = gsum$species[order(gsum$weight)]) # sort by weight

# Horizontal bars

h1 <- ggplot(gsum, aes(x = species, y = bs, fill = nest))+
  geom_bar(stat = "identity", width=0.85, colour="black")+
  coord_flip()+
  geom_errorbar(aes(ymin = gsum$ymin, ymax = gsum$ymax),width=0.85/2, position=position_dodge(.9))+
  scale_y_continuous(expand = c(0, 0), limits = c(0,0.95))+
  scale_x_discrete(breaks=unique(gsum$species), 
                   labels=gsum$species) +
  scale_fill_manual(values=c("#bebebe", "#FFFFFF"),
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
  geom_point(x = 10, y = 0.37, color = "black", shape=23, fill="black", size=5) +
  geom_text(aes(label = paste0(weight, "g"), y = 0.9, hjust = 1, family = "serif"), size = 5)


#ggsave("plot_NS_3.png", h1, height = 7, width = 7, dpi = 600)

# Mixed effects models for Gough species
# Explore possible outliers
a <- gDat %>%
  filter(nest == "b") # Below-ground nesters; repeat for above-ground and each season
quantiles <- quantile(a$bs, probs = c(.25, .75))
range <- 1.5 * IQR(a$bs)
range

norm <- subset(a, a$bs > (quantiles[1] - range)
               & a$bs < (quantiles[2] + range))

# Rescale weights to facilitate model convergence
linMap <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from}
gDat$wscale <- linMap(gDat$weight, 0, 1)

# transform breeding success so that model residuals are normal and variances are equal
gDat$bs2 <- sqrt(gDat$bs) 
df <- data.frame(species = unique(gDat$species), s2 = 1:10)
gDat <- merge(gDat, df, by = "species")

wr.fit1 <- glmer(bs2 ~ (1|y2), family = binomial, data = gDat)
wr.fit2 <- glmer(bs2 ~ s2 + (1|y2), family = binomial, data = gDat)
wr.fit3 <- glmer(bs2 ~ nest +  (1|y2), family = binomial, data = gDat)
wr.fit4 <- glmer(bs2 ~ season + (1|y2), family = binomial, data = gDat)
wr.fit5 <- glmer(bs2 ~ wscale + (1|y2), family = binomial, data = gDat)
wr.fit6 <- glmer(bs2 ~ s2 + nest + (1|y2), family = binomial, data = gDat)
wr.fit7 <- glmer(bs2 ~ s2 + season + (1|y2), family = binomial, data = gDat)
wr.fit8 <- glmer(bs2 ~ s2 + wscale + (1|y2), family = binomial, data = gDat)
wr.fit9 <- glmer(bs2 ~ nest + season + (1|y2), family = binomial, data = gDat)
wr.fit10 <- glmer(bs2 ~ nest + wscale + (1|y2), family = binomial, data = gDat)
wr.fit11 <- glmer(bs2 ~ season + wscale + (1|y2), family = binomial, data = gDat)
wr.fit12 <- glmer(bs2 ~ s2 + nest + season + (1|y2), family = binomial, data = gDat)
wr.fit13 <- glmer(bs2 ~ s2 + nest + wscale + (1|y2), family = binomial, data = gDat)
wr.fit14 <- glmer(bs2 ~ s2 + season + wscale + (1|y2), family = binomial, data = gDat)
wr.fit15 <- glmer(bs2 ~ nest + season + wscale + (1|y2), family = binomial, data = gDat)
wr.fit16 <- glmer(bs2 ~ s2 + nest + season + wscale + (1|y2), family = binomial, data = gDat, na.action=na.fail)
anova(wr.fit5, wr.fit8)
histogram(residuals(wr.fit5))
qqnorm(residuals(wr.fit5))

# Use broom to extract model components and save into data frame objects
# Extract variable stats using tidy, e.g.:
# wr.fit5 %>% tidy
# And/or model data using:
# wr.fit5 %>% augment() %>% head(3)
# Include likelihood ratio tests (lr)
# Model names
modcomp <- data.frame(mod = sprintf("wr.fit%d",seq(1:16)),
                      rbind(glance(wr.fit1),glance(wr.fit2),glance(wr.fit3),
                            glance(wr.fit4),glance(wr.fit5),glance(wr.fit6),
                            glance(wr.fit7),glance(wr.fit8),glance(wr.fit9),
                            glance(wr.fit10),glance(wr.fit11),glance(wr.fit12),
                            glance(wr.fit13),glance(wr.fit14),glance(wr.fit15),
                            glance(wr.fit16)),
                      lr = c((AIC(wr.fit1)-AIC(wr.fit1))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit2))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit3))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit4))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit5))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit6))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit7))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit8))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit9))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit10))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit11))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit12))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit13))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit14))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit15))*log2(exp(1)),
                             (AIC(wr.fit1)-AIC(wr.fit16))*log2(exp(1))
                      ))
# Compare models
modcomp

# Dredge all data
wr.aic <- dredge(wr.fit16, rank = "AIC")
var.imp <- importance(wr.aic) # Variable  -importance weights
top.set <- subset(wr.aic, delta <2) # Top set of models (<2 delta AIC)
top.mod <- subset(wr.aic, delta == 0) # Best apporoximatng model

# Average model
wr.avg <- model.avg(top.set)
summary(wr.avg)
confint(mbDat.avg, level = 0.95)



# Mixed effects model between Gough-analogue pairs
anDat$bs2 <- sqrt(anDat$bs) # tranformation for normality

ga.fit1 <- glmer(bs ~ (1|y2), family=binomial, data = anDat)
ga.fit2 <- glmer(bs ~ l2 + (1|y2), family=binomial, data = anDat)
anova(ga.fit1, ga.fit2)




# Comparative productivity
# Read summed & paired Gough-analogue csv
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
names(gDAT) <- c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.gough.bs"="gough.bs")
gDAT["min_err"] <- minmax(dat = gDAT, m = gDAT$gough.bs, s = bDAT$gough.se, l=0, fn="min.e")
gDAT["max_err"] <- minmax(dat = gDAT, m = gDAT$gough.bs, s = bDAT$gough.se, u = 1, fn="max.e")
gDAT["obs.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, pe = "obs.c")
gDAT["min.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, min.e = gDAT$min_err, pe = "min.c")
gDAT["max.c"] <- popfunc (dat = gDAT, pop = gDAT$pop, bs = gDAT$gough.bs, max.e = gDAT$max_err, pe = "max.c")
gDAT[9,c(6:8)] <- c(150, 110, 190) # Southern giant petrel; rounded to 50, a-priori
gDAT[10,c(6:8)] <- c(550, 400, 700) # Tristan albatross; rounded to 50, a-priori
gDAT["mean.c"] <- ((gDAT$max.c-gDAT$obs.c)+
                     (gDAT$obs.c-gDAT$min.c))/2


# Typical data (uninvaded islands)

tDAT <- data.frame(bDAT$species, bDAT$est_pop, bDAT$typical.bs)
names(tDAT) <- c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.typical.bs"="typical.bs")

tDAT["min_err"] <- minmax(dat = tDAT, m = tDAT$typical.bs, s = bDAT$typical.se, l=0, fn="min.e")
tDAT["max_err"] <- minmax(dat = tDAT, m = tDAT$typical.bs, s = bDAT$typical.se, u=1, fn="max.e")
tDAT["obs.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, pe = "obs.c")
tDAT["min.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, min.e = tDAT$min_err, pe = "min.c")
tDAT["max.c"] <- popfunc (dat = tDAT, pop = tDAT$pop, bs = tDAT$typical.bs, max.e = tDAT$max_err, pe = "max.c")
tDAT[9,c(7:9)] <- c(180, 150, 200) # Southern giant petrel; rounded to 50, a-priori
tDAT[10,c(7:9)] <- c(1250, 1150, 1400) # Tristan albatross; rounded to 50, a-priori
tDAT["mean.c"] <- ((tDAT$max.c-tDAT$obs.c)+
                     (tDAT$obs.c-tDAT$min.c))/2


# Calculate difference between population estimates and relative impact of mice

pDAT <- data.frame(bDAT$species, bDAT$est_pop, bDAT$bs.diff)
names(pDAT) <- c("bDAT.species"="species", "bDAT.est_pop" = "pop",
                       "bDAT.bs.diff"="difference.bs")
pDAT["min_err"] <- minmax(dat = pDAT, m = pDAT$difference.bs, s = bDAT$se.diff, l=0, fn="min.e")
pDAT["max_err"] <- minmax(dat = pDAT, m = pDAT$difference.bs, s = bDAT$se.diff, u=1, fn="max.e")
pDAT["obs.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, pe = "obs.c")
pDAT["min.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, min.e = pDAT$min_err, pe = "min.c")
pDAT["max.c"] <- popfunc (dat = pDAT, pop = pDAT$pop, bs = pDAT$difference.bs, max.e = pDAT$max_err, pe = "max.c")
pDAT[9,c(7:9)] <- c(30, 0, 50) # Southern giant petrel; rounded to 50, a-priori
pDAT[10,c(7:9)] <- c(700, 650, 750) # Tristan albatross; rounded to 50, a-priori
pDAT["mean.c"] <- ((pDAT$max.c-pDAT$obs.c)+
                     (pDAT$obs.c-pDAT$min.c))/2

# Write to csv

write.csv(gDAT, file = "d_gough.csv")
write.csv(tDAT, file = "d_typical.csv")
write.csv(pDAT, file = "d_population.csv")

# Perform Chi-square test

pop <- data.frame(gDAT$obs.c, tDAT$obs.c) 
chisq.test(pop)
