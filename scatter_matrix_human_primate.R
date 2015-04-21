library(ggplot2)
library(reshape2)
library(colorspace)
library(scales)
library(gridExtra)

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

get_plot <- function(x_lab, y_lab, dataset, min_cp, max_cp, colors, pchs) {
	p1 <- ggplot(dataset, aes_string(x=x_lab, y=y_lab)) + geom_point(aes(colour = code, pch = code), size=1.5, alpha = 0.7) + 
			  theme_bw() + theme(panel.grid.minor = element_blank()) + 
			  theme(legend.position = "none", axis.title = element_blank()) +
			  scale_colour_manual(name = title, values=colors) + scale_shape_manual(name = title, values=pchs) +
			  scale_x_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp)) + 
			  scale_y_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp))
	return(p1)
}

get_plot_with_legend <- function(x_lab, y_lab, dataset, min_cp, max_cp, colors, pchs) {
    p1 <- ggplot(dataset, aes_string(x=x_lab, y=y_lab)) + geom_point(aes(colour = code, pch = code), size=5, alpha = 0.7) +
              theme_bw() + theme(panel.grid.minor = element_blank()) +
              scale_colour_manual(name = title, values=colors) + scale_shape_manual(name = title, values=pchs) +
              scale_x_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp)) +
              scale_y_continuous(limits=c(floor(min_cp),max_cp), breaks=seq(floor(min_cp),max_cp))
    return(p1)
}


theme_bare <- theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(), 
  axis.ticks.margin = unit(c(0,0,0,0), "lines"), 
  legend.position = "none", 
  panel.background = element_rect(fill = "white"), 
  panel.border = element_blank(), 
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(), 
  panel.margin = unit(c(0,0,0,0), "lines"), 
  plot.background = element_rect(fill = "white"),
  plot.margin = unit(c(0,0,0,0), "lines")
)

args <- commandArgs(trailingOnly = T)
infile <- args[1]
outfile <- args[2]
dataset <- args[3]

regions <- args[4:length(args)]

if(length(regions) < 2) {
	print("At least two regions required")
	quit()
}

dat = read.table(infile, sep="\t", header=T)
dat <- dat[dat$name %in% regions,]
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
nhumans = length(humans) - 2
colors.human <- c(grDevices::rainbow(nhumans), "gray20", "black")
colors.nhp <- grDevices::gray(seq(0, 0.75, length.out=4))
pchs.human <- c(rep(1, nhumans), 3, 4)
pchs.nhp <- rep(0, 4)

dat.sample$species <- sapply(dat.sample$code, function(x) if(x %in% humans) return("human") else return("non human primate"))

dat.len = dim(dat.sample)[1]

dat.combined <- dat.sample[with(dat.sample, order(code)),]

#dat.combined$pch <- sapply(dat.combined$code, function(x) as.factor(pchs[match(x, unique(dat.combined$code))]))
#dat.combined$color <- sapply(dat.combined$code, function(x) as.factor(colors.final[match(x, unique(dat.combined$code))]))

min_cp = min(dat.combined[,regions])
#max_cp = min(8, max(dat.combined[,regions]))
max_cp = max(dat.combined[,regions])

p1_with_legend <- get_plot_with_legend(regions[1], regions[2], dat.combined, min_cp, max_cp, c(colors.human, colors.nhp), c(pchs.human, pchs.nhp))
plot.legend <- g_legend(p1_with_legend)

plot.list = list()

for(i in 1:length(regions)) {
	for(j in 1:length(regions)) {
		if(i == j) {
			p1 <- ggplot(dat.combined) + annotate("text", x=(max_cp - min_cp)/2, y=(max_cp-min_cp)/2, label=regions[i], size=10) + theme_bw() + theme_bare
		} else if(i > j) {
			p1 <- get_plot(regions[j], regions[i], dat.combined[dat.combined$species == "human",], min_cp, max_cp, colors.human, pchs.human)
		} else {
			p1 <- get_plot(regions[j], regions[i], dat.combined[dat.combined$species == "non human primate",], min_cp, max_cp, colors.nhp, pchs.nhp)
		}
		plot.list[[(i-1)*length(regions) + j]] <- p1
	}
}

width = max(6, 5*length(regions) + 2)

plots <- do.call(arrangeGrob, c(plot.list, ncol=length(regions)))
g <- arrangeGrob(plots, plot.legend, ncol=2, widths=c(0.65*width, 0.25*width))
print(g)
ggsave(g, filename=outfile, height=max(6, 5*length(regions)), width=width)
