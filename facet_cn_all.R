library(ggplot2)
library(reshape2)
library(colorspace)
library(scales)
library(plyr)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
dataset <- args[3]

comparisons <- args[4:length(args)]
dat = read.table(infile, sep="\t", header=T)

coords = c(as.vector(unique(dat$field)), "\ncode")
title = paste(coords, collapse="\n")

dat.sample <- dcast(dat, sample + super_pop + code ~ name, value.var = "copy_num")

### Settings for HGDP or 1kg, both with Archaics and NHP 
if(dataset == "hgdp") {
	dat.sample$code <- factor(dat.sample$code, levels=c("AFR_African", "AMR_Americas", "OCN_Oceanic", "EA_East_Asian", "SA_South_Asian", "SIB_Siberian", "WEA_Western_European", "ARC_Archaics", "D&N_Denisova_and_Neanderthal", "CHIMPANZEE", "BONOBO", "GORILLA", "ORANGUTAN"))
	humans <- c("AFR_African", "AMR_Americas", "OCN_Oceanic", "EA_East_Asian", "SA_South_Asian", "SIB_Siberian", "WEA_Western_European", "ARC_Archaics", "D&N_Denisova_and_Neanderthal")
} else if(dataset == "1kg") {
	dat.sample$code <- factor(dat.sample$code, levels=c("AFR_African", "AMR_Admixed_American", "EAS_East_Asian", "EUR_European", "SAS_South_Asian", "ARC_Archaics", "D&N_Denisova_and_Neanderthal", "CHIMPANZEE", "BONOBO", "GORILLA", "ORANGUTAN"))
	humans <- c("AFR_African", "AMR_Admixed_American", "EAS_East_Asian", "EUR_European", "SAS_South_Asian", "ARC_Archaics", "D&N_Denisova_and_Neanderthal")
} else {
	print(paste("dataset", dataset, "not recognized. Exiting...", sep=" "))
	quit(1)
}
###	Get colors and symbols
ncolors = length(humans) - 2
colors.final <- c(grDevices::rainbow(ncolors), "gray20", "black", grDevices::gray(seq(0, 0.75, length.out=4)))
pchs <- c(rep(1, ncolors), 3, 4, rep(0, 4))

dat.sample$species <- sapply(dat.sample$code, function(x) if(x %in% humans) return("human") else return("non human primate"))

dat.list = list()
dat.len = dim(dat.sample)[1]

for(i in 1:(length(comparisons)/2)) {
	dat.tmp <- dat.sample
	dat.tmp$first <- dat.sample[,comparisons[2*i-1]]
	dat.tmp$second <- dat.sample[,comparisons[2*i]]
	dat.tmp$type <- paste(comparisons[2*i-1], "vs", comparisons[2*i], sep=" ")
	dat.list[[i]] <- dat.tmp
}

dat.combined <- ldply(dat.list, data.frame)
dat.combined <- dat.combined[with(dat.combined, order(code)),]

dat.combined$pch <- sapply(dat.combined$code, function(x) as.factor(pchs[match(x, unique(dat.combined$code))]))

min_cp = min(c(dat.combined$first, dat.combined$second))
max_cp = min(max(c(dat.combined$first, dat.combined$second)))

p_combined <- ggplot(dat.combined, aes(first, second)) + geom_point(aes(colour = code, pch = code), size=1.5, alpha = 0.7) + 
			  xlab("ARHGAP11A copy number") + ylab("ARHGAP11B copy number") + theme_bw() + theme(panel.grid.minor = element_blank()) + 
			  facet_grid(type ~ species) + guides(color=guide_legend(override.aes = list(size=5, pch=pchs), title=title)) +
			  scale_colour_manual(name = title, values=colors.final) + scale_shape_manual(name = title, values=pchs) +
			  scale_x_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp)) + 
			  scale_y_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp))
			  #scale_fill_manual(values=alpha(colors.final, alphas))

ggsave(p_combined, filename=outfile, height=max(4, length(comparisons)), width=8)
