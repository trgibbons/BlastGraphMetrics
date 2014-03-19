#!/usr/bin/Rscript --vanilla --default-packages=utils

library(tools)
library(grDevices)
library(ggplot2)

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
                     levels=c("BitScore","BitPerResidue","BitScoreRatio",
                              "Evalue", "p(Evalue)"))
kpc$Metric <- factor(kpc$Metric,
                     levels=c("BitScore","BitPerResidue","BitScoreRatio",
                              "Evalue", "p(Evalue)"))

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
                 aes(x=Inflation, y=ClustersPerKOG))
cpk.gg <- cpk.gg +theme_bw()
cpk.gg <- cpk.gg +facet_grid(.~Metric)
cpk.gg <- cpk.gg +coord_trans(y="log10")
cpk.gg <- cpk.gg +ggtitle(bquote(atop("Clusters per KOG", atop(.(args[1]), ""))))
cpk.gg <- cpk.gg +ylab("log10(Clusters per KOG)")
cpk.gg <- cpk.gg +xlab("MCL Inflation Parameter")

                
# Plot heatmaps to PDFs
pdf(paste(cpk.pref,"_I11-50_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(cpk.gg +geom_tile(color="grey42", aes(fill=cpk$ClusterCount))
             +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()

pdf(paste(cpk.pref,"_I11-50_doublelog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(cpk.gg +geom_tile(color="grey42", aes(fill=log10(cpk$ClusterCount)))
             +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()


# Create ggplot object for KOGs per Cluster statistics
kpc.gg <- ggplot(data=kpc, aes(x=Inflation, y=KOGsPerCluster))
kpc.gg <- kpc.gg +geom_tile(color="grey42")
kpc.gg <- kpc.gg +theme_bw()
kpc.gg <- kpc.gg +facet_grid(.~Metric)
kpc.gg <- kpc.gg +coord_trans(y="log10")
kpc.gg <- kpc.gg +ggtitle(bquote(atop("KOGs per Cluster", atop(.(args[2]), ""))))
kpc.gg <- kpc.gg +ylab("log10(KOGs per Clusters)")
kpc.gg <- kpc.gg +xlab("MCL Inflation Parameter")

# Plot heatmaps to PDFs
pdf(paste(kpc.pref,"_I11-50_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(kpc.gg +geom_tile(color="grey42", aes(fill=kpc$ClusterCount))
             +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()

pdf(paste(kpc.pref,"_I11-50_doublelog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(kpc.gg +geom_tile(color="grey42", aes(fill=log10(kpc$ClusterCount)))
             +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()


# Subset data to focus on just inflation parameter values between [1.1,2]
cpk12 <- droplevels.data.frame(subset(cpk, as.numeric(Inflation)<=2.0, drop=TRUE))
kpc12 <- droplevels.data.frame(subset(kpc, as.numeric(Inflation)<=2.0, drop=TRUE))

# ...then repeat the whole thing
# Create ggplot object for Clusters per KOG statistics
cpk12.gg <- ggplot(data=cpk12, aes(x=Inflation, y=ClustersPerKOG))
cpk12.gg <- cpk12.gg +theme_bw()
cpk12.gg <- cpk12.gg +facet_grid(.~Metric)
cpk12.gg <- cpk12.gg +coord_trans(y="log10")
cpk12.gg <- cpk12.gg +ggtitle(bquote(atop("Clusters per KOG", atop(.(args[1]), ""))))
cpk12.gg <- cpk12.gg +ylab("log10(Clusters per KOG)")
cpk12.gg <- cpk12.gg +xlab("MCL Inflation Parameter")

# Plot heatmaps to PDFs
pdf(paste(cpk.pref,"_I11-20_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(cpk12.gg +geom_tile(color="grey42", aes(fill=cpk12$ClusterCount))
               +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()

pdf(paste(cpk.pref,"_I11-20_doublelog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(cpk12.gg +geom_tile(color="grey42", aes(fill=log10(cpk12$ClusterCount)))
               +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()


# Create ggplot object for KOGs per Cluster statistics
kpc12.gg <- ggplot(data=kpc12,
                   aes(x=Inflation, y=KOGsPerCluster))
kpc12.gg <- kpc12.gg +theme_bw()
kpc12.gg <- kpc12.gg +facet_grid(.~Metric)
kpc12.gg <- kpc12.gg +coord_trans(y="log10")
kpc12.gg <- kpc12.gg +ggtitle(bquote(atop("KOGs per Cluster", atop(.(args[2]), ""))))
kpc12.gg <- kpc12.gg +ylab("log10(KOGS per Cluster)")
kpc12.gg <- kpc12.gg +xlab("MCL Inflation Parameter")

# Plot heatmaps to PDFs
pdf(paste(kpc.pref,"_I11-20_semilog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(kpc12.gg +geom_tile(color="grey42", aes(fill=kpc12$ClusterCount))
               +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()

pdf(paste(kpc.pref,"_I11-20_doublelog_8.5x3.pdf", sep=""), width=8.5, height=3)
print(kpc12.gg +geom_tile(color="grey42", aes(fill=log10(kpc12$ClusterCount)))
               +scale_fill_gradientn(colours=rp, name="Cluster Count"))
dev.off()




