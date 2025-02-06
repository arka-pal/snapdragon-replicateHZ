library(data.table)
library(tidyverse)

#### read data ----
AvePlaFRYe_w10KBs10KB = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('./AvePla-FRYe/AveFR-AveY-PlaFR-PlaY_Chr',chrom,'_byPos_w10000_s10000.csv.gz')
  dat = fread(file)
  return(dat)
}
str(AvePlaFRYe_w10KBs10KB)
dat = AvePlaFRYe_w10KBs10KB[sites>=20]


#### Quantiles ----
q99_AveFR_AveY = quantile(dat$Fst_AveFR_AveY, probs = c(0.5,0.95,0.99))[3]
q95_AveFR_AveY = quantile(dat$Fst_AveFR_AveY, probs = c(0.5,0.95,0.99))[2]
q99_PlaFR_PlaY = quantile(dat$Fst_PlaFR_PlaY, probs = c(0.5,0.95,0.99))[3]
q95_PlaFR_PlaY = quantile(dat$Fst_PlaFR_PlaY, probs = c(0.5,0.95,0.99))[2]
q99_AveFR_PlaFR = quantile(dat$Fst_AveFR_PlaFR, probs = c(0.5,0.95,0.99))[3]
q95_AveFR_PlaFR = quantile(dat$Fst_AveFR_PlaFR, probs = c(0.5,0.95,0.99))[2]
q99_AveY_PlaY = quantile(dat$Fst_AveY_PlaY, probs = c(0.5,0.95,0.99))[3]
q95_AveY_PlaY = quantile(dat$Fst_AveY_PlaY, probs = c(0.5,0.95,0.99))[2]
q99_AveFR_PlaY = quantile(dat$Fst_AveFR_PlaY, probs = c(0.5,0.95,0.99))[3]
q95_AveFR_PlaY = quantile(dat$Fst_AveFR_PlaY, probs = c(0.5,0.95,0.99))[2]
q99_AveY_PlaFR = quantile(dat$Fst_AveY_PlaFR, probs = c(0.5,0.95,0.99))[3]
q95_AveY_PlaFR = quantile(dat$Fst_AveY_PlaFR, probs = c(0.5,0.95,0.99))[2]


#### Define colors ----
dat = dat %>% 
  mutate(windowCol_col = case_when(
    Fst_AveFR_AveY >= q95_AveFR_AveY & Fst_PlaFR_PlaY >= q95_PlaFR_PlaY ~ "black",
    Fst_AveFR_AveY >= q95_AveFR_AveY ~ "grey70",
    Fst_PlaFR_PlaY >= q95_PlaFR_PlaY ~ "grey70",
    TRUE ~ "grey90"),
    windowCol_geo = case_when(
      Fst_AveFR_PlaFR >= q95_AveFR_PlaFR & Fst_AveY_PlaY >= q95_AveY_PlaY ~ "black",
      Fst_AveFR_PlaFR >= q95_AveFR_PlaFR ~ "grey70",
      Fst_AveY_PlaY >= q95_AveY_PlaY ~ "grey70",
      TRUE ~ "grey90"),
    windowCol_colgeo = case_when(
      Fst_AveFR_PlaY >= q95_AveFR_PlaY & Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "black",
      Fst_AveFR_PlaY >= q95_AveFR_PlaY ~ "grey70",
      Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "grey70",
      TRUE ~ "grey90"),
    
    ## Correlations with AveY
    windowCol_2_6 = case_when(
      Fst_AveFR_AveY >= q95_AveFR_AveY & Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "black",
      Fst_AveFR_AveY >= q95_AveFR_AveY ~ "grey70",
      Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "grey70",
      TRUE ~ "grey90"),
    
    windowCol_4_6 = case_when(
      Fst_AveY_PlaY >= q95_AveY_PlaY & Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "black",
      Fst_AveY_PlaY >= q95_AveY_PlaY ~ "grey70",
      Fst_AveY_PlaFR >= q95_AveY_PlaFR ~ "grey70",
      TRUE ~ "grey90"),
    
    windowCol_2_4 = case_when(
      Fst_AveFR_AveY >= q95_AveFR_AveY & Fst_AveY_PlaY >= q95_AveY_PlaY ~ "black",
      Fst_AveFR_AveY >= q95_AveFR_AveY ~ "grey70",
      Fst_AveY_PlaY >= q95_AveY_PlaY ~ "grey70",
      TRUE ~ "grey90")
  )

## YELLOW
cremosa = c('Chr1', 821293, 1095277)
aurina = c('Chr2', 951685, 1153653)
flavia = c('Chr2', 53368062, 53569435)
sulf = c('Chr4', 38248572, 38450336)

## MAGENTA
rubia = c('Chr5', 6169966, 6372151)
rosel = c('Chr6', 52789323, 53162505)


#### plots ----


## ------
par(xaxs='i', yaxs='i', mfrow=c(1,3))

## within zones
plot(Fst_AveFR_AveY~Fst_PlaFR_PlaY, dat, las = 1, main='', xlab='', ylab='', pch=20, cex=1, xlim=c(0,0.5), ylim=c(0,1), col=dat$windowCol_col)

alphaVal = 1; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 25
points(Fst_AveFR_AveY~Fst_PlaFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(h=c(q95_AveFR_AveY, q99_AveFR_AveY), col=c('black', 'red'), lty=1:2)
abline(v=c(q95_PlaFR_PlaY, q99_PlaFR_PlaY), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')


## ------
## between zones
plot(Fst_AveY_PlaY~Fst_AveFR_PlaFR, dat, las = 1, main='',
     xlab='', ylab='', pch=20, cex=1,
     xlim=c(0,0.5), ylim=c(0,1), col=dat$windowCol_geo)

alphaVal = 0.8; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 25
points(Fst_AveY_PlaY~Fst_AveFR_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(h=c(q95_AveY_PlaY, q99_AveY_PlaY), col=c('black', 'red'), lty=1:2)
abline(v=c(q95_AveFR_PlaFR, q99_AveFR_PlaFR), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')


## ------ ##
## between zones and vars
plot(Fst_AveY_PlaFR~Fst_AveFR_PlaY, dat, las = 1, main='',
     xlab='', ylab='', pch=20, cex=1,
     xlim=c(0,0.5), ylim=c(0,1), col=dat$windowCol_colgeo)

alphaVal = 0.8; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 24
points(Fst_AveY_PlaFR~Fst_AveFR_PlaY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(h=c(q95_AveY_PlaFR, q99_AveY_PlaFR), col=c('black', 'red'), lty=1:2)
abline(v=c(q95_AveFR_PlaY, q99_AveFR_PlaY), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')




## plots AveY ----
par(xaxs='i', yaxs='i', mfrow=c(1,3))
## 2-6 comparison
plot(Fst_AveY_PlaFR~Fst_AveFR_AveY, dat, las = 1, #main='Color comparison',
     xlab='', ylab='', pch=20, cex=1,
     xlim=c(0,1), ylim=c(0,1), col=dat$windowCol_2_6)

alphaVal = 0.8; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 25
points(Fst_AveY_PlaFR~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(v=c(q95_AveFR_AveY, q99_AveFR_AveY), col=c('black', 'red'), lty=1:2)
abline(h=c(q95_AveY_PlaFR, q99_AveY_PlaFR), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')


###
## 2-4 comparison
plot(Fst_AveY_PlaY~Fst_AveFR_AveY, dat, las = 1, #main='Geography comparison',
     xlab='', ylab='', pch=20, cex=1,
     xlim=c(0,1), ylim=c(0,1), col=dat$windowCol_2_4)

alphaVal = 0.8; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 25
points(Fst_AveY_PlaY~Fst_AveFR_AveY, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(v=c(q95_AveFR_AveY, q99_AveFR_AveY), col=c('black', 'red'), lty=1:2)
abline(h=c(q95_AveY_PlaY, q99_AveY_PlaY), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')


###
## 6-4 comparison
plot(Fst_AveY_PlaY~Fst_AveY_PlaFR, dat, las = 1, #main='Geography comparison',
     xlab='', ylab='', pch=20, cex=1,
     xlim=c(0,1), ylim=c(0,1), col=dat$windowCol_4_6)

alphaVal = 0.8; geneSize = 1.2
geneName = rubia; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = cremosa; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = aurina; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)
geneName = sulf; geneColor = 'blue'; geneShape = 20
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       col=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = rosel; geneColor = 'magenta'; geneShape = 24
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

alphaVal = 1; geneSize = 1.8
geneName = flavia; geneColor = 'yellow'; geneShape = 25
points(Fst_AveY_PlaY~Fst_AveY_PlaFR, 
       dat[scaffold==geneName[1] & start>=as.integer(geneName[2]) & end<=as.integer(geneName[3])], 
       bg=alpha(geneColor,alphaVal), pch=geneShape, cex=geneSize)

abline(v=c(q95_AveY_PlaFR, q99_AveY_PlaFR), col=c('black', 'red'), lty=1:2)
abline(h=c(q95_AveY_PlaY, q99_AveY_PlaY), col=c('black', 'red'), lty=1:2)
abline(a=0, b=1, lty='longdash')

