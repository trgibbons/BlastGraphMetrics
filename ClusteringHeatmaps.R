#!/usr/bin/Rscript --vanilla --default-packages=utils

library(tools)
library(grDevices)
library(ggplot2)

# Get command line arguments
args <- commandArgs(TRUE)

# Color palette
rp <- c("lightblue","yellow","orange","red","purple")

# Read in data
cpk <- read.table(args[1], header=TRUE)
kpc <- read.table(args[2], header=TRUE)


# Create ggplot object for Clusters per KOG statistics
cpk$Metric <- factor(cpk$Metric,
                     levels=c("BitScore","BitPerResidue","BitScoreRatio",
                              "Evalue"))
cpk.gg <- ggplot(data=cpk,
                 aes(x=Inflation, y=ClustersPerKOG, fill=KOGCount))
cpk.gg <- cpk.gg +geom_tile(color="grey42")
cpk.gg <- cpk.gg +theme_bw()
cpk.gg <- cpk.gg +facet_grid(.~Metric)
cpk.gg <- cpk.gg +theme(axis.title.x=element_blank(),
                        axis.title.y=element_blank())
cpk.gg <- cpk.gg +theme(strip.text.x=element_blank(),
                        strip.background=element_blank())
cpk.gg <- cpk.gg +scale_fill_gradientn(name=element_blank(), colours=rp,
                                       limits=c(1,458))


# Create ggplot object for KOGs per Cluster statistics
kpc$Metric <- factor(kpc$Metric,
                     levels=c("BitScore","BitPerResidue","BitScoreRatio",
                              "Evalue"))
kpc.gg <- ggplot(data=kpc,
                 aes(x=Inflation, y=KOGsPerCluster,fill=ClusterCount))
kpc.gg <- kpc.gg +geom_tile(color="grey42")
kpc.gg <- kpc.gg +theme_bw()
kpc.gg <- kpc.gg +facet_grid(.~Metric)
kpc.gg <- kpc.gg +theme(axis.title.x=element_blank(),
                        axis.title.y=element_blank())
kpc.gg <- kpc.gg +theme(strip.text.x=element_blank(),
                        strip.background=element_blank())
kpc.gg <- kpc.gg +scale_fill_gradientn(name=element_blank(), colours=rp,
                                       limits=c(1,600))


# Get prefixes for output files
cpk.pref <- file_path_sans_ext(args[1])
kpc.pref <- file_path_sans_ext(args[2])


# Plot heatmaps
pdf(paste(cpk.pref,"_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(cpk.gg +coord_trans(y="log10"))
dev.off()

pdf(paste(kpc.pref,"_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(kpc.gg +coord_trans(y="log10"))
dev.off()

