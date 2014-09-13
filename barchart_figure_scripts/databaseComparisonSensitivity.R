#!/usr/bin/Rscript --vanilla --default-packages=utils

library(tools)
library(grDevices)
library(plyr)
library(ggplot2)
library(RColorBrewer)

# BlastGraphMetrics directory
bgmdir <- "/Volumes/OdinsSaddlebag/Research/2014/clustering/BlastGraphMetrics"

# Read in data
cegma.fn <- "cegma_1/shf_rnd/1e-5/nrm_dmnd/cegma_111111_shf_rnd_1e-5_nrm_dmnd_clusters_per_kog_summary.Rtab"
cegma <- read.table(paste(bgmdir, cegma.fn, sep="/"), header=TRUE)
cegma$DataSet <- as.factor("CEGMA")
cegma$ClusterCount <- cegma$ClusterCount/100

eck.fn <- "eck_11111111111/shf_rnd/1e-5/nrm_dmnd/eck_11111111111_shf_rnd_1e-5_nrm_dmnd_clusters_per_kog_summary.Rtab"
eck <- read.table(paste(bgmdir, eck.fn, sep="/"), header=TRUE)
eck$DataSet <- as.factor("ECK")
eck$ClusterCount <- eck$ClusterCount/100

kog.fn <- "kog_1/shf_rnd/1e-5/nrm_dmnd/kog_1111111_shf_rnd_1e-5_nrm_dmnd_clusters_per_kog_summary.Rtab"
kog <- read.table(paste(bgmdir, kog.fn, sep="/"), header=TRUE)
kog$DataSet <- as.factor("KOG")
kog$ClusterCount <- kog$ClusterCount/1000

dsdf <- rbind(cegma, eck, kog)


# Refactor the Metric columns to set the order
dsdf$Metric <- factor(dsdf$Metric,
                      levels=c("-Log10Evalue", "BitScore", "BitScoreRatio",
                               "AnchoredLength"))

dsdf <- ddply(dsdf, c("Metric", "Inflation"),
              transform, TotesClusters=sum(ClusterCount))

dsdf$Legend <- as.character(dsdf$ClustersPerKOG)
dsdf$Legend[as.numeric(dsdf$Legend) >= 5] <- "5+"
dsdf$Legend <- as.factor(dsdf$Legend)

# If you're lucky enough to not be red/green color blind, I think this color scheme looks better
# right <- "#006d2c"
# wrong <- brewer.pal(name="YlOrRd", n=7)[4:7]
right <- "#4575b4"
wrong <- rev(brewer.pal(name="RdYlBu", n=11))[7:10]

dsdf.reds <- nlevels(dsdf$Legend)-1
if (dsdf.reds > 0) {
  dsdf.clrs = c(right, wrong[1:nlevels(dsdf$Legend)-1])
} else {
  dsdf.clrs = c(right)
}

# Create ggplot object for Clusters per KOG statistics
dsdf.gg <- ggplot(data=dsdf,
                  aes(x=Inflation, weight=ClusterCount))
dsdf.gg <- dsdf.gg +geom_bar(aes(fill=Legend), binwidth=0.1)
dsdf.gg <- dsdf.gg +scale_fill_manual(
                        values=dsdf.clrs,
                        guide=guide_legend(
                            title="Clusters\nper KOG",
                            title.hjust=0.5))
dsdf.gg <- dsdf.gg +facet_grid(DataSet~Metric, scales="free", space="free")
dsdf.gg <- dsdf.gg +theme_bw()
dsdf.gg <- dsdf.gg +ggtitle("Sensitivity")
dsdf.gg <- dsdf.gg +ylab("KOG Count")
dsdf.gg <- dsdf.gg +xlab("MCL Inflation Parameter")
dsdf.gg <- dsdf.gg +geom_hline(
                        yintercept=seq(1.0,max(dsdf$ClusterCount),1.0),
                        color="lightgray", size=0.05)

# Plot barcharts to PDFs
pdf("datasetComparisonSensitivity.pdf", width=8.5, height=6)
dsdf.gg
dev.off()
