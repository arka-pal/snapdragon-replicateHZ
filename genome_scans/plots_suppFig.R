library(data.table)
library(tidyverse)
library(foreach)
library(gridExtra)
# library(patchwork)
library(viridis)
setwd('~/snap_hap_repHZ/genome_scans/')
source('~/snap_hap_repHZ/genome_scans/_scripts/functions_genomescans.R')
plot_theme <- theme_bw() + theme(panel.grid = element_blank())


## Read popgen params ----
AvePlaFRYe_w10KBs10KB = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('./AvePla-FRYe/AveFR-AveY-PlaFR-PlaY_Chr',chrom,'_byPos_w10000_s10000.csv.gz')
  dat = fread(file)
  return(dat)
}
AvePlaFRYe_w10KBs10KB = modify_for_ManhattanPlot(AvePlaFRYe_w10KBs10KB)
axisToPlot = axis_for_ManhattanPlot(AvePlaFRYe_w10KBs10KB, 'grey10','grey80')
str(AvePlaFRYe_w10KBs10KB)


## Read recRate ----
rec_AveM = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('~/snap_hap_repHZ/recRate_bySnps/Chr',chrom,'/AveM/chr',chrom,'_AveM_recRate_w10Ks10K.txt')
  dat = fread(file)
  return(dat)
}
rec_AveY = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('~/snap_hap_repHZ/recRate_bySnps/Chr',chrom,'/AveY/chr',chrom,'_AveY_recRate_w10Ks10K.txt')
  dat = fread(file)
  return(dat)
}
rec_PlaM = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('~/snap_hap_repHZ/recRate_bySnps/Chr',chrom,'/PlaM/chr',chrom,'_PlaM_recRate_w10Ks10K.txt')
  dat = fread(file)
  return(dat)
}
rec_PlaY = foreach(chrom=c(1:8), .combine=rbind) %do% {
  file=paste0('~/snap_hap_repHZ/recRate_bySnps/Chr',chrom,'/PlaY/chr',chrom,'_PlaY_recRate_w10Ks10K.txt')
  dat = fread(file)
  return(dat)
}

# Add chrom_start_end coloum
AvePlaFRYe_w10KBs10KB = AvePlaFRYe_w10KBs10KB %>% mutate(chrom_start_end = paste(scaffold, start, end, sep='_')) 
commonCols = c('chrom', 'windowStart', 'windowEnd', 'windowMid', 'sites')
recRate = cbind(rec_PlaY[,..commonCols], rec_PlaM[,3], rec_PlaY[,3], rec_AveM[,3], rec_AveY[,3])
colnames(recRate) = c('chrom', 'windowStart', 'windowEnd', 'windowMid', 'rec_sites', 'rec_PlaM', 'rec_PlaY', 'rec_AveM', 'rec_AveY')
recRate = recRate %>% mutate(chrom_start_end = paste(chrom, windowStart, windowEnd, sep='_')) 
str(recRate)
str(AvePlaFRYe_w10KBs10KB)

# Combine the 2 datasets
combinedDat = merge(recRate, AvePlaFRYe_w10KBs10KB, by='chrom_start_end')
str(combinedDat)

# Subset windows with >20 sites
dat = combinedDat[sites >= 20]
summary(dat)


## Color palletes ----
viridisPallete = viridis(100, option = 'viridis', direction = -1)
magmaPallete = viridis(100, option = 'magma', direction = -1)
cividisPallete = viridis(100, option = 'cividis', direction = +1)
pointColors = function(dataToColor, maxValue, colorPallete){
  scale = maxValue/100
  z_scaled = dataToColor %/% scale
  z_scaled = pmax(pmin(z_scaled, 100), 1)
  return(colorPallete[z_scaled])
}

par(mfrow=c(1,3))
color_min = 0; color_max = 1
color_legend = seq(color_min, color_max, length.out = 100)

image(1, color_legend, t(matrix(color_legend, ncol = 1)), 
      col = viridis(100, option = 'viridis', direction = -1), axes=F, 
      xaxt = "n", yaxt = "n", xlab='', ylab='')
axis(2, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)
axis(4, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)

image(1, color_legend, t(matrix(color_legend, ncol = 1)), 
      col = viridis(100, option = 'magma', direction = -1), axes=F, 
      xaxt = "n", yaxt = "n", xlab='', ylab='')
axis(2, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)
axis(4, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)

image(1, color_legend, t(matrix(color_legend, ncol = 1)), 
      col = viridis(100, option = 'cividis', direction = 1), axes=F, 
      xaxt = "n", yaxt = "n", xlab='', ylab='')
axis(2, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)
axis(4, at=seq(0,1,0.5), lab=NA, tcl=1, col=NA, col.ticks = 'white', las = 1)

## Biplots of pi, Dxy, Fst, rec ----
par(mfrow=c(3,2))

##Pi/Pi, col: Fst
xLim=c(0,0.5); yLim=c(0,0.5); point_shape = 1; alphaVal = 1; maxVal = 0.7 #fst

datSorted = dat %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "Pi within PlaM", ylab = 'Pi within PlaY')
abline(a=0, b=1, lty='longdash')

datSorted = dat %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal),
     xlab = "Pi within AveM", ylab = 'Pi within AveY')
abline(a=0, b=1, lty='longdash')


##meanPi/Fst, col: Dxy
xLim=c(0,0.5); yLim=c(0,0.8); point_shape = 1; alphaVal = 1; maxVal = 0.66 #dxy

datSorted = dat %>% arrange(dxy_PlaFR_PlaY)
plot((datSorted$pi_PlaFR + datSorted$pi_PlaY)/2, datSorted$Fst_PlaFR_PlaY,
     xlim=xLim, ylim=yLim, pch = point_shape, xlab = "Mean Pi", ylab = 'Fst',
     col=alpha(pointColors(datSorted$dxy_PlaFR_PlaY, maxVal, magmaPallete), alphaVal))

datSorted = dat %>% arrange(dxy_AveFR_AveY)
plot((datSorted$pi_AveFR + datSorted$pi_AveY)/2, datSorted$Fst_AveFR_AveY,
     xlim=xLim, ylim=yLim, pch = point_shape, xlab = "Mean Pi", ylab = 'Fst',
     col=alpha(pointColors(datSorted$dxy_AveFR_AveY, maxVal, magmaPallete), alphaVal))

## meanRec/Fst, col: meanPi
xLim=c(0,175); yLim=c(0,0.8); point_shape = 1; alphaVal = 1; maxVal = 0.4 #pi

datSorted = dat %>% mutate(meanPi = rowMeans(select(., pi_PlaFR, pi_PlaY))) %>% arrange(meanPi) %>% select(-meanPi)
plot((datSorted$rec_PlaM + datSorted$rec_PlaY)/2, datSorted$Fst_PlaFR_PlaY,
     xlim=xLim, ylim=yLim, pch = point_shape, xlab = "Mean Recombination Rate", ylab = 'Fst',
     col=alpha(pointColors((datSorted$pi_PlaFR + datSorted$pi_PlaY)/2, maxVal, cividisPallete), alphaVal))

datSorted = dat %>% mutate(meanPi = rowMeans(select(., pi_AveFR, pi_AveY))) %>% arrange(meanPi) %>% select(-meanPi)
plot((datSorted$rec_AveM + datSorted$rec_AveY)/2, datSorted$Fst_AveFR_AveY,
     xlim=xLim, ylim=yLim, pch = point_shape, xlab = "Mean Recombination Rate", ylab = 'Fst',
     col=alpha(pointColors((datSorted$pi_AveFR + datSorted$pi_AveY)/2, maxVal, cividisPallete), alphaVal))



## Biplots of pi per chromosome (Planoles) ----
par(mfrow=c(2,4), mar=c(4,4,2,2), cex=1, xaxs='i', yaxs='i', las=1)
xLim=c(0,0.5); yLim=c(0,0.5); point_shape = 20; alphaVal = 1; maxVal = 0.7 #fst

#Chr1
datSorted = dat[chrom == 'Chr1'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr2
datSorted = dat[chrom == 'Chr2'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr3
datSorted = dat[chrom == 'Chr3'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr4
datSorted = dat[chrom == 'Chr4'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr5
datSorted = dat[chrom == 'Chr5'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr6
datSorted = dat[chrom == 'Chr6'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr7
datSorted = dat[chrom == 'Chr7'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr8
datSorted = dat[chrom == 'Chr8'] %>% arrange(Fst_PlaFR_PlaY)
plot(pi_PlaY~pi_PlaFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_PlaFR_PlaY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')



## Biplots of pi per chromosome (Avellanet) ----
par(mfrow=c(2,4), mar=c(4,4,2,2), cex=1, xaxs='i', yaxs='i', las=1)
xLim=c(0,0.5); yLim=c(0,0.5); point_shape = 20; alphaVal = 1; maxVal = 0.7 #fst

#Chr1
datSorted = dat[chrom == 'Chr1'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr2
datSorted = dat[chrom == 'Chr2'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr3
datSorted = dat[chrom == 'Chr3'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr4
datSorted = dat[chrom == 'Chr4'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr5
datSorted = dat[chrom == 'Chr5'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr6
datSorted = dat[chrom == 'Chr6'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr7
datSorted = dat[chrom == 'Chr7'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')

#Chr8
datSorted = dat[chrom == 'Chr8'] %>% arrange(Fst_AveFR_AveY)
plot(pi_AveY~pi_AveFR, datSorted, xlim=xLim, ylim=yLim, pch = point_shape, 
     col=alpha(pointColors(datSorted$Fst_AveFR_AveY, maxVal, viridisPallete), alphaVal), 
     xlab = "", ylab = '')
abline(a=0, b=1, lty='longdash')


## Blowup Chr1 & 2 ----
par(mfrow=c(2,4), mar=c(2,5,2,2), mgp=c(2.5,1,0), cex.lab=1.2, yaxs='i')

dat_chr1 = dat[chrom == 'Chr1'] %>% arrange(mid)
span = 0.02
fst_Ave_smooth = loess.smooth.weights(dat_chr1$mid, dat_chr1$Fst_AveFR_AveY, sites=dat_chr1$sites, span=span, family='gaussian')
dxy_Ave_smooth = loess.smooth.weights(dat_chr1$mid, dat_chr1$dxy_AveFR_AveY, sites=dat_chr1$sites, span=span, family='gaussian')
pi_AveM_smooth = loess.smooth.weights(dat_chr1$mid, dat_chr1$pi_AveFR, sites=dat_chr1$sites, span=span, family='gaussian')
pi_AveY_smooth = loess.smooth.weights(dat_chr1$mid, dat_chr1$pi_AveY, sites=dat_chr1$sites, span=span, family='gaussian')
rec_Ave_smooth = loess.smooth.weights(dat_chr1$mid, (dat_chr1$rec_AveM+dat_chr1$rec_AveY)/2, sites=dat_chr1$sites, span=span, family='gaussian')
#Fst
plot(Fst_AveFR_AveY~mid, dat_chr1, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(F)[ST]), ylim=c(0,0.8))
lines(fst_Ave_smooth~mid, dat_chr1, cex=0.8, col='black', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,1,0.4), cex.axis=1)
#Dxy
plot(dxy_AveFR_AveY~mid, dat_chr1, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(D)[XY]), ylim=c(0,0.8))
lines(dxy_Ave_smooth~mid, dat_chr1, col='blue', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,1,0.4), cex.axis=1)
#Pi AveM
plot(pi_AveFR~mid, dat_chr1, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(pi)[w]), ylim=c(0,0.5))
lines(pi_AveM_smooth~mid, dat_chr1, col='magenta', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,0.5,0.2), cex.axis=1)
#Pi AveY
plot(pi_AveY~mid, dat_chr1, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(pi)[w]), ylim=c(0,0.5))
lines(pi_AveY_smooth~mid, dat_chr1, col='gold', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,0.5,0.2), cex.axis=1)


dat_chr2 = dat[chrom == 'Chr2'] %>% arrange(mid)
span = 0.02
fst_Ave_smooth = loess.smooth.weights(dat_chr2$mid, dat_chr2$Fst_AveFR_AveY, sites=dat_chr2$sites, span=span, family='gaussian')
dxy_Ave_smooth = loess.smooth.weights(dat_chr2$mid, dat_chr2$dxy_AveFR_AveY, sites=dat_chr2$sites, span=span, family='gaussian')
pi_AveM_smooth = loess.smooth.weights(dat_chr2$mid, dat_chr2$pi_AveFR, sites=dat_chr2$sites, span=span, family='gaussian')
pi_AveY_smooth = loess.smooth.weights(dat_chr2$mid, dat_chr2$pi_AveY, sites=dat_chr2$sites, span=span, family='gaussian')
rec_Ave_smooth = loess.smooth.weights(dat_chr2$mid, (dat_chr2$rec_AveM+dat_chr2$rec_AveY)/2, sites=dat_chr2$sites, span=span, family='gaussian')

#Fst
plot(Fst_AveFR_AveY~mid, dat_chr2, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(F)[ST]), ylim=c(0,0.8))
lines(fst_Ave_smooth~mid, dat_chr2, cex=0.8, col='black', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,1,0.4), cex.axis=1)
#Dxy
plot(dxy_AveFR_AveY~mid, dat_chr2, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(D)[XY]), ylim=c(0,0.8))
lines(dxy_Ave_smooth~mid, dat_chr2, col='blue', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,1,0.4), cex.axis=1)
#Pi AveM
plot(pi_AveFR~mid, dat_chr2, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(pi)[w]), ylim=c(0,0.5))
lines(pi_AveM_smooth~mid, dat_chr2, col='magenta', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,0.5,0.2), cex.axis=1)

#Pi AveY
plot(pi_AveY~mid, dat_chr2, cex=0.8, col=alpha('grey', 0.5), 
     xaxt='n', yaxt='n', type='p', pch=20, xlab='', ylab=expression(italic(pi)[w]), ylim=c(0,0.5))
lines(pi_AveY_smooth~mid, dat_chr2, col='gold', lwd=3)
axis(1, at=seq(0,70e6,20e6), lab=seq(0,70,20), cex.axis=1)
axis(2, at=seq(0,0.5,0.2), cex.axis=1)

