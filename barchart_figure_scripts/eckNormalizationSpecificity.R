#!/usr/bin/Rscript --vanilla --default-packages=utils

library(tools)
library(grDevices)
library(plyr)
library(ggplot2)
library(RColorBrewer)

# BlastGraphMetrics directory
bgmdir <- "/Volumes/OdinsSaddlebag/Research/2014/clustering/BlastGraphMetrics"

# Data sets
dss <- c("raw","nrm")

# Read in data
for (ds in dss) {
  fn <- paste("eck_11111111111/shf_rnd/1e-5/",ds,"_dmnd/eck_11111111111_shf_rnd_1e-5_",ds,"_dmnd_kogs_per_cluster_summary.Rtab", sep="")
  tdf <- read.table(paste(bgmdir, fn, sep="/"), header=TRUE)
  tdf$DataSet <- as.factor(ds)
  if (exists("eck")) {
    eck <- rbind(eck, tdf)
  } else {
    eck <- tdf
  }
}

# Refactor the Metric columns to set the order
eck$Metric <- factor(eck$Metric,
                     levels=c("-Log10Evalue", "BitScore", "BitScoreRatio",
                              "AnchoredLength"))

eck <- ddply(eck, c("Metric", "Inflation"),
             transform, TotesClusters=sum(ClusterCount))

eck$Legend <- as.character(eck$KOGsPerCluster)
eck$Legend[as.numeric(eck$Legend) >= 5] <- "5+"
eck$Legend <- as.factor(eck$Legend)

# If you're lucky enough to not be red/green color blind, I think this color scheme looks better
# right <- "#006d2c"
# wrong <- brewer.pal(name="YlOrRd", n=7)[4:7]
right <- "#4575b4"
wrong <- rev(brewer.pal(name="RdYlBu", n=11))[7:10]

eck.reds <- nlevels(eck$Legend)-1
if (eck.reds > 0) {
  eck.clrs = c(right, wrong[1:nlevels(eck$Legend)-1])
} else {
  eck.clrs = c(right)
}

# Create ggplot object for ECKs per Cluster statistics
eck.gg <- ggplot(data=eck,
                 aes(x=Inflation, weight=ClusterCount))
eck.gg <- eck.gg +geom_bar(aes(fill=Legend), binwidth=0.1)
eck.gg <- eck.gg +scale_fill_manual(
                      values=eck.clrs,
                      guide=guide_legend(
                          title="ECKs per\nCluster",
                          title.hjust=0.5))
eck.gg <- eck.gg +facet_grid(Normalization~Metric, scales="free", space="free")
eck.gg <- eck.gg +theme_bw()
eck.gg <- eck.gg +ggtitle("Specificity")
eck.gg <- eck.gg +ylab("Cluster Count")
eck.gg <- eck.gg +xlab("MCL Inflation Parameter")
eck.gg <- eck.gg +geom_hline(yintercept=seq(458,max(458,max(eck$ClusterCount)),458),
                             linetype="dashed")
eck.gg <- eck.gg +geom_hline(yintercept=seq(100,max(eck$ClusterCount),100),
                             color="lightgray", size=0.05)

# Plot barcharts to PDFs
pdf("eckNormalizationSpecificity.pdf", width=8.5, height=4)
eck.gg
dev.off()
