#!/usr/bin/Rscript --vanilla --default-packages=utils

library(tools)
library(grDevices)
library(plyr)
library(ggplot2)
library(RColorBrewer)

# Get command line arguments
args <- commandArgs(TRUE)

# Get prefixes for output files
cpk.pref <- file_path_sans_ext(args[1])
kpc.pref <- file_path_sans_ext(args[2])

# Color palette
rp <- c("yellow","green","purple","red")

# Read in data
cpk <- read.table(args[1], header=TRUE)
kpc <- read.table(args[2], header=TRUE)

# Refactor the Metric columns to set the order
cpk$Metric <- factor(cpk$Metric,
                     levels=c("-Log10Evalue", "BitScore", "BitScoreRatio",
                              "AnchoredLength"))
kpc$Metric <- factor(kpc$Metric,
                     levels=c("-Log10Evalue", "BitScore", "BitScoreRatio",
                              "AnchoredLength"))

cpk <- ddply(cpk, c("Metric", "Inflation"),
             transform, TotesClusters=sum(ClusterCount))
kpc <- ddply(kpc, c("Metric", "Inflation"),
             transform, TotesClusters=sum(ClusterCount))

cpk$Legend <- as.character(cpk$ClustersPerKOG)
cpk$Legend[as.numeric(cpk$Legend) >= 5] <- "5+"
cpk$Legend <- as.factor(cpk$Legend)

kpc$Legend <- as.character(kpc$KOGsPerCluster)
kpc$Legend[as.numeric(kpc$Legend) >= 5] <- "5+"
kpc$Legend <- as.factor(kpc$Legend)

# If you're lucky enough to not be red/green color blind, I think this color scheme looks better
# right <- "#006d2c"
# wrong <- brewer.pal(name="YlOrRd", n=7)[4:7]
right <- "#4575b4"
wrong <- rev(brewer.pal(name="RdYlBu", n=11))[7:10]

cpk.reds <- nlevels(cpk$Legend)-1
if (cpk.reds > 0) {
    cpk.clrs = c(right, wrong[1:nlevels(cpk$Legend)-1])
} else {
    cpk.clrs = c(right)
}

kpc.reds <- nlevels(cpk$Legend)-1
if (kpc.reds > 0) {
    kpc.clrs = c(right, wrong[1:nlevels(kpc$Legend)-1])
} else {
    kpc.clrs = c(right)
}


# For future reference:
#  To remove labels before exporting to Illustrator
# gg <- gg +theme(axis.title.x=element_blank(),
#                 axis.title.y=element_blank())
# gg <- gg +theme(strip.text.x=element_blank(),
#                 strip.background=element_blank())
#
# To choose numbers to be indicated on legends
# cpk.brks.lin <- c(1,200,400,458,500,max(cpk$ClusterCount))
# cpk.brks.log <- c(log10(1),log10(200),log10(400),log10(458),log10(500),
#                  log10(max(cpk$ClusterCount)))
#                                   
# kpc.brks.lin <- c(1,200,400,458,500,max(kpc$ClusterCount))
# kpc.brks.log <- c(log10(1),log10(200),log10(400),log10(458),log10(500),
#                  log10(max(kpc$ClusterCount)))


# Create ggplot object for Clusters per KOG statistics
cpk.gg <- ggplot(data=cpk,
                 aes(x=Inflation, weight=ClusterCount))
cpk.gg <- cpk.gg +geom_bar(aes(fill=Legend), binwidth=0.1)
cpk.gg <- cpk.gg +scale_fill_manual(
                      values=cpk.clrs,
                      guide=guide_legend(
                                title="Clusters\nper ECK",
                                title.hjust=0.5))
cpk.gg <- cpk.gg +facet_grid(.~Metric)
cpk.gg <- cpk.gg +theme_bw()
cpk.gg <- cpk.gg +ggtitle(bquote(atop("Sensitivity",
                                      atop(.(args[1]), ""))))
cpk.gg <- cpk.gg +ylab("KOG Count")
cpk.gg <- cpk.gg +xlab("MCL Inflation Parameter")
cpk.gg <- cpk.gg +geom_hline(yintercept=seq(458,max(458,max(cpk$ClusterCount)),458),
                             linetype="dashed")
cpk.gg <- cpk.gg +geom_hline(yintercept=seq(100,max(cpk$ClusterCount),100),
                             color="lightgray", size=0.05)
# cpk.gg <- cpk.gg +geom_vline(xintercept=seq(1,4),
#                              color="lightgray", size=0.05)
                
# Plot barcharts to PDFs
pdf(paste(cpk.pref,"_I1-6_barcharts_8.5x3.pdf", sep=""), width=8.5, height=3)
cpk.gg
dev.off()

pdf(paste(cpk.pref,"_I1-3_barcharts_8.5x3.pdf", sep=""), width=8.5, height=3)
cpk.gg +xlim(1, 3)
dev.off()


# Create ggplot object for Clusters per KOG statistics
kpc.gg <- ggplot(data=kpc,
                 aes(x=Inflation, weight=ClusterCount))
kpc.gg <- kpc.gg +geom_bar(aes(fill=Legend), binwidth=0.1)
kpc.gg <- kpc.gg +scale_fill_manual(
                      values=kpc.clrs,
                      guide=guide_legend(
                                title="KOGs per\nCluster",
                                title.hjust=0.5))
kpc.gg <- kpc.gg +facet_grid(.~Metric)
kpc.gg <- kpc.gg +theme_bw()
kpc.gg <- kpc.gg +ggtitle(bquote(atop("Specificity",
                                      atop(.(args[2]), ""))))
kpc.gg <- kpc.gg +ylab("Cluster Count")
kpc.gg <- kpc.gg +xlab("MCL Inflation Parameter")
kpc.gg <- kpc.gg +geom_hline(yintercept=seq(458,max(458,max(kpc$ClusterCount)),458),
                             linetype="dashed")
kpc.gg <- kpc.gg +geom_hline(yintercept=seq(100,max(kpc$ClusterCount),100),
                             color="lightgray", size=0.1)
# kpc.gg <- kpc.gg +geom_vline(xintercept=seq(1,4),
#                              color="lightgray", size=0.1)

# Plot heatmaps to PDFs
pdf(paste(kpc.pref,"_I1-6_barcharts_8.5x3.pdf", sep=""), width=8.5, height=3)
kpc.gg
dev.off()

pdf(paste(kpc.pref,"_I1-3_barcharts_8.5x3.pdf", sep=""), width=8.5, height=3)
kpc.gg +xlim(1, 3)
dev.off()
