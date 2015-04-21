library(ggplot2)
library(gridExtra)

g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

args <- commandArgs(trailingOnly=TRUE)
input.raw <- args[1]
input.GMM <- args[2]
output.prefix <- args[3]
region.name <- args[4]
output.type <- args[5]
plot.title <- args[6]

copy_nums <- read.table(input.raw, header=TRUE, sep="\t")
copy_nums <- copy_nums[copy_nums$name==region.name,]

gmm.copy_nums <- read.table(input.GMM, header=TRUE, sep="\t")
gmm.copy_nums <- gmm.copy_nums[gmm.copy_nums$name==region.name,]
gmm.copy_nums$code <- factor(gmm.copy_nums$code, levels = rev(unique(gmm.copy_nums$code)), ordered = TRUE)
sorted.gmm.copy_nums <-gmm.copy_nums[with(gmm.copy_nums, order(code, pop)),]
sorted.gmm.copy_nums$pop <- factor(sorted.gmm.copy_nums$pop, levels=unique(sorted.gmm.copy_nums$pop))
sorted.gmm.copy_nums$super_pop <- factor(sorted.gmm.copy_nums$super_pop, levels=unique(sorted.gmm.copy_nums$super_pop), ordered = TRUE)

copy_nums$code <- factor(copy_nums$code, levels = rev(unique(copy_nums$code)), ordered = TRUE)
legend_cn = copy_nums

sorted.copy_nums <- copy_nums[with(copy_nums, order(code, pop)),]
sorted.copy_nums$pop <- factor(sorted.copy_nums$pop, levels=unique(sorted.copy_nums$pop))
sorted.copy_nums$super_pop <- factor(sorted.copy_nums$super_pop, levels=unique(sorted.copy_nums$super_pop), ordered = TRUE)

sd.count <- sd(sorted.copy_nums$copy_num)
min.count <- floor(min(sorted.copy_nums$copy_num))
max.count <- ceiling(max(sorted.copy_nums$copy_num))

threshold <- mean(copy_nums$copy_num) + 8*sd(copy_nums$copy_num)

xlab = "Super population"
ylab = "Copy number"

plot.legend <- ggplot(sorted.gmm.copy_nums, aes(x=code, y=copy_num, fill=code, name="Super population")) + 
        geom_violin() + theme_bw() + xlab(xlab) + ylab(ylab) + 
        scale_y_continuous(breaks=0:max.count, limits=c(-0.5, max.count+0.5)) + 
		guides(fill = guide_legend(reverse = TRUE))

pop.legend <- g_legend(plot.legend)

p2 <- ggplot(sorted.gmm.copy_nums, aes(factor(copy_num), fill=super_pop)) + geom_bar() + theme_bw() + xlab(ylab) + ylab("Count") + theme(legend.position="none")

p1 <- ggplot(sorted.copy_nums, aes(x=super_pop, y=copy_num, fill=code)) + 
geom_point(alpha=0.5, colour='black', solid=T, size=1, position = position_jitter(h = 0, w=0.1)) + 
theme_bw() + coord_flip() + xlab(xlab) + ylab(ylab) + theme(legend.position="none") + 
scale_y_continuous(breaks=0:max.count, limits=c(-0.5, max.count+0.5), minor_breaks=c())

if(output.type == "pdf") {
    pdf(output.prefix, width=12, height=3)
    grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
    dev.off()
} else if(output.type == "png") {
png(output.prefix, width=800, height=200*3)
grid.arrange(arrangeGrob(p1, p2, ncol=2), pop.legend, ncol=2, widths=c(4/5, 1/5), main=plot.title)
dev.off()
} else print(paste("Unsupported file type", output.type))
