while(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
while(!require(ggplot2)) install.packages("ggplot2", repos = "http://cran.us.r-project.org")
while(!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org")


##################################################################################

#plotting parameters

#plot height and width
height              <- 40
width               <- 40
#ratio of tree to poplabels
ratio               <- c(10,4)

#tree linewidth
tree_lwd            <- 3
#mutation size
mut_size            <- 10
#size of | indicating population label
poplabels_shapesize <- 10 
#population label text size
poplabels_textsize  <- 30

##################################################################################

### Function: filetype()
filetype <- function(path){
    f = file(path)
    ext = summary(f)$class
    close.connection(f)
    ext
}

### Function: Treeview()
TreeView <- function(filename_plot, years_per_gen, axis_textsize, ...){

  plotcoords      <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
  plotcoords[3:4] <- plotcoords[3:4] * years_per_gen

  p <- ggplot()          + geom_segment(data = subset(plotcoords, seg_type != "m"), aes(x = x_begin, xend = x_end, y = y_begin, yend = y_end), ...) +
                            #geom_point(data = subset(plotcoords, seg_type == "m"), aes(x = x_begin, y = y_begin), color = "black", size = 2) +
                            theme(text = element_text(size=axis_textsize),
                                  axis.line=element_blank(),axis.text.x=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.title.x=element_blank(),
                                  axis.title.y=element_blank(),
                                  # legend.position="top",
                                  legend.position="none",
                                  panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
                                  panel.grid.minor=element_blank(), plot.background=element_blank(), 
                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"), 
                                  # legend.key = element_blank(),
                                  # legend.key.width= unit(3, "line"),
                                  # legend.key.height= unit(1.5, "line"),
                                  # legend.text=element_text(size=35),
                                  strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 30, unit = "pt")) + 
                             scale_x_continuous(limits = c(0, max(plotcoords$x_begin)+1)) + scale_y_continuous(limits = c(0, max(plotcoords$y_end)))

  return(p)

}

### Function: AddMutations()
AddMutations <- function(filename_plot, filename_mut, years_per_gen, ...){

	plotcoords      <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
	plotcoords[3:4] <- plotcoords[3:4] * years_per_gen

	mut_on_branches     <- read.table(paste(filename_plot,".plotcoords.mut", sep = ""), header = T)

	muts <- subset(plotcoords, seg_type == "v" | seg_type == "t")
	muts <- merge(mut_on_branches, muts, by = "branchID")

	muts %>% group_by(branchID) %>% mutate(y_begin = (1:length(y_begin)) * (max(y_end) - min(y_begin))/(1+length(y_begin)) + min(y_begin), y_end = y_begin ) -> muts
	mut <- read.table(filename_mut, skip = 1, sep = ";")
	mut <- mut[,c(2,8)]
	colnames(mut) <- c("pos", "is_flipped")
	mut$is_flipped <- as.factor(mut$is_flipped)
	muts <- merge(muts, mut, by = "pos")
	p <- geom_point(data = muts, aes(x = x_begin, y = y_begin, colour = is_flipped), ...)

	return(p)

}


### Function: PopLabels()
# PopLabels <- function(filename_plot, filename_poplabels, text_size = 100, ...){

#   plotcoords <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
#   poplabels  <- read.table(filename_poplabels, header = T)[,2:4]

#   tips <- subset(plotcoords, seg_type == "t")
# 	if(all(is.na(poplabels[,3])) || any(poplabels[,3] != 1)){
# 		tips <- cbind(tips, population = poplabels[ceiling((tips$branchID+1)/2),1], region = poplabels[ceiling((tips$branchID+1)/2),2])
# 	}else{
# 		tips <- cbind(tips, population = poplabels[tips$branchID+1,1], region = poplabels[tips$branchID+1,2])
# 	}
#   unique_region <- unique(poplabels[,2])
#   p <- ggplot() + geom_point(data = tips, aes(x = x_begin, y = population, color = population), ...) +
#                    theme(text = element_text(size=text_size),
#                                   axis.line=element_blank(),axis.text.x=element_blank(),
#                                   axis.ticks=element_blank(),
#                                   axis.title.x=element_blank(),
#                                   axis.title.y=element_blank(),legend.position="bottom",
#                                   panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#                                   panel.grid.minor=element_blank(),plot.background=element_blank(), 
#                                   strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                                   strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 60, unit = "pt")) +
#                    guides(color = F, shape = F) + 
#                    scale_x_continuous(limits = c(0, max(tips$x_begin)+1))

#   return(p)

# }
PopLabels <- function(filename_plot, filename_poplabels, 
                      shape, poplabels_shapesize, poplabels_textsize, ...){
    plotcoords <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
    poplabels  <- read.table(filename_poplabels, header = T)[,2:4]

    tips <- subset(plotcoords, seg_type == "t")
    if(all(is.na(poplabels[,3])) || any(poplabels[,3] != 1)){
        tips <- cbind(tips, population = poplabels[ceiling((tips$branchID+1)/2),1], region = poplabels[ceiling((tips$branchID+1)/2),2])
    } else {
        tips <- cbind(tips, population = poplabels[tips$branchID+1,1], region = poplabels[tips$branchID+1,2])
    }
    unique_region <- unique(poplabels[,2])

    p <-    tips %>% 
            # arrange(region, x_begin, population) %>% 
            # mutate(population = fct_reorder(population, x_begin)) %>%
            ggplot(aes(x = x_begin, y = population, color = population)) + 
                geom_point(pch=shape, size=poplabels_shapesize) + 
                theme(text = element_text(face='italic', size=poplabels_textsize),
                                  axis.line=element_blank(),axis.text.x=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.title.x=element_blank(),
                                  axis.title.y=element_blank(),legend.position="bottom",
                                  panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                                  panel.grid.minor=element_blank(),plot.background=element_blank(), 
                                  strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                                  strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 60, unit = "pt")) +
                guides(color = 'none', shape = 'none') + 
                scale_color_manual(values = rep(c('magenta', 'gold3'),2)) +
                scale_x_continuous(limits = c(0, max(tips$x_begin)+1))
    
    return(p)
}
#################################################################################

### Function: drawtipLabels()
drawTipLabels = function(filename_plot, filename_poplabels, tipInfo){
    plotcoords <- read.table(paste(filename_plot,".plotcoords", sep = ""), header = T)
    poplabels  <- read.table(filename_poplabels, header = T)[,2:4]
    tips <- subset(plotcoords, seg_type == "t")
    if(all(is.na(poplabels[,3])) || any(poplabels[,3] != 1)){
        tips <- cbind(tips, population = poplabels[ceiling((tips$branchID+1)/2),1], region = poplabels[ceiling((tips$branchID+1)/2),2])
    } else {
        tips <- cbind(tips, population = poplabels[tips$branchID+1,1], region = poplabels[tips$branchID+1,2])
    }
    tips$color = tipInfo[match(as.integer(tips$branchID), tipInfo$IDs), "tipColors"]
    
    tipLabels = ggplot(tips) + 
                    geom_text(aes(x = x_begin, y = 0, label = branchID, color=color), angle=90, size=3) +
                    labs(y = 'Tips') +
                    scale_color_identity() +
                    theme(text = element_text(face='italic', size=poplabels_textsize),
                          axis.line=element_blank(), 
                          axis.text.x=element_blank(),
                          # axis.text.y=element_blank(),
                          axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          # axis.title.y=element_blank(),legend.position="bottom",
                          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),plot.background=element_blank(), 
                          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                          strip.text = element_text(face="bold"), plot.margin = margin(t = 0, r = 20, b = 60, l = 60, unit = "pt"))

    return(tipLabels)
}


#################################################################################
drawTree = function(PATH_TO_RELATE='~/_softwares/relate_v1.2.2', 
                    filename_haps, filename_sample, filename_anc, filename_mut,
                    filename_poplabels='./AvePla.MY.n74.poplabels',
                    years_per_gen, snp, filename_plot,
                    geneName, makeFiles = FALSE, 
                    tree_lwd=0.5, axis_textsize=15, mut_size=5, 
                    shape='|', poplabels_shapesize=8, poplabels_textsize=20, 
                    ratio=c(10,4)){
    
    if(makeFiles == TRUE){
        ## run RelateTreeView to extract tmp files for plotting tree
        system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeView --anc ", filename_anc, " --mut ", filename_mut, " --snp_of_interest ", as.integer(snp), " -o ", filename_plot))
        system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_anc, " --mut ", filename_mut, " --haps ", filename_haps, " --sample ", filename_sample, " --snp_of_interest ", as.integer(snp), " -o ", filename_plot))
    }


    ## Plot tree
    p1 <- TreeView(filename_plot, years_per_gen, lwd = tree_lwd, axis_textsize) + 
            AddMutations(filename_plot, filename_mut, years_per_gen, size = mut_size) #+ 
            #scale_y_continuous(trans = "log10")
    
    # # some modifications to theme
    # p1 <- p1 + theme(axis.text.y = element_text(size = rel(2.3)), 
    #                  legend.title = element_text(size = rel(1)), 
    #                  legend.text = element_text(size = rel(1)))  +  
    #            scale_color_manual(labels = c("unflipped", "flipped"), values = c("red", "blue"), drop = FALSE) +
    #            guides(color = guide_legend(nrow = 2, title = "")) 

    ## plot population label
    p2 <- PopLabels(filename_plot, filename_poplabels, 
                    shape, poplabels_shapesize, poplabels_textsize)

    finalPlot = plot_grid(p1, p2, rel_heights = ratio, labels = "", align = "v", ncol = 1)
    return(finalPlot)
}

#################################################################################
drawTreeWithTips = function(PATH_TO_RELATE='~/_softwares/relate_v1.2.2', 
                    filename_haps, filename_sample, filename_anc, filename_mut,
                    filename_poplabels='./AvePla.MY.n74.poplabels',
                    years_per_gen, snp, filename_plot,
                    geneName, tipInfo, makeFiles = FALSE, 
                    tree_lwd=0.5, axis_textsize=15, mut_size=5, 
                    shape='|', poplabels_shapesize=8, poplabels_textsize=20, 
                    ratio=c(10,2,4)){
    
    if(makeFiles == TRUE){
        ## run RelateTreeView to extract tmp files for plotting tree
        system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode TreeView --anc ", filename_anc, " --mut ", filename_mut, " --snp_of_interest ", as.integer(snp), " -o ", filename_plot))
        system(paste0(PATH_TO_RELATE, "/bin/RelateTreeView --mode MutationsOnBranches --anc ", filename_anc, " --mut ", filename_mut, " --haps ", filename_haps, " --sample ", filename_sample, " --snp_of_interest ", as.integer(snp), " -o ", filename_plot))
    }


    ## Plot tree
    p1 <- TreeView(filename_plot, years_per_gen, lwd = tree_lwd, axis_textsize) + 
            AddMutations(filename_plot, filename_mut, years_per_gen, size = mut_size) #+ 
            #scale_y_continuous(trans = "log10")
    
    # # some modifications to theme
    # p1 <- p1 + theme(axis.text.y = element_text(size = rel(2.3)), 
    #                  legend.title = element_text(size = rel(1)), 
    #                  legend.text = element_text(size = rel(1)))  +  
    #            scale_color_manual(labels = c("unflipped", "flipped"), values = c("red", "blue"), drop = FALSE) +
    #            guides(color = guide_legend(nrow = 2, title = "")) 

    p2 = drawTipLabels(filename_plot, filename_poplabels, tipInfo)
    
    ## plot population label
    p3 <- PopLabels(filename_plot, filename_poplabels, 
                    shape, poplabels_shapesize, poplabels_textsize)

    finalPlot = plot_grid(p1, p2, p3, rel_heights = ratio, labels = "", align = "v", ncol = 1)
    return(finalPlot)
}


##################################################################################

# pdf(paste(filename_plot,".pdf", sep = ""), height = height, width = width)
# plot_grid(p1,p2, rel_heights = ratio, labels = "", align = "v", ncol = 1)
# dev.off()

# # delete tmp files for plotting tree
# system(paste("rm ",filename_plot, ".plotcoords", sep = ""))
# system(paste("rm ",filename_plot, ".plotcoords.mut", sep = ""))