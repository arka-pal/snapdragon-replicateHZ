library(data.table)
library(tidyverse)
library(foreach)

## -----
modify_for_ManhattanPlot = function(dat){
    dat %>% 
        # Compute chromosome size
        group_by(scaffold) %>% summarise(chromLength = max(end)) %>% 
    
        # Calculate cumulative position of each chromosome
        mutate(total = cumsum(chromLength)-chromLength) %>% select(-chromLength) %>%
    
        # Add this info to the initial dataset
        left_join(dat, ., by=c("scaffold"="scaffold")) %>%
    
        # Add a cumulative position of each SNP
        arrange(scaffold, end) %>%
        mutate(midCum = mid+total)
}

## -----
axis_for_ManhattanPlot = function(modified_dat, col1, col2){
    modified_dat %>%
        group_by(scaffold) %>% summarize(chromCenter = (max(midCum, na.rm=T) + min(midCum, na.rm=T))/2) %>%
        mutate(col = rep(c(col1, col2), 4))
}

## -----
loess.smooth.weights = function(position, value, sites, span, family){
    y.loess <- loess(value ~ position, span = span, weights = sites, family=family)
    y.predict <- predict(y.loess, position)
    y.predict
}