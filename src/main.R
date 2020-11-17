# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import libraries
library(ggplot2)
library(pheatmap)
library(vegan)
library(stringr)
library(dplyr)

# import metadata
metadata_ibd <- read.csv("data/metadata_ibd.txt", header = T, stringsAsFactors = F)

# import data
mpa.non_sorting <- read.table("data/mpa.non_sorting.txt", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
mpa.org <- mpa.non_sorting[which(!is.na(str_extract(colnames(mpa.non_sorting), "s__")))]
colnames(mpa.org) <- sub("s__", "", sapply(str_split(colnames(mpa.org), "\\."), function(x) x[7]))

# filtering data
Numzz <- colSums(mpa.org==0)
Savezz <-  (Numzz<nrow(mpa.org)*0.75)
mpa.org.sbs <- mpa.org[Savezz]
mpa.org.sbs <- mpa.org.sbs[order(colSums(mpa.org.sbs), decreasing = T)]

ids <- c(metadata_ibd$df_sample, metadata_ibd$donor_id_manual)
ids <- unique(ids[!is.na(ids)])

mpa.org.sbs <- mpa.org.sbs[rownames(mpa.org.sbs) %in% ids,]
mpa.org.sbs <- mpa.org.sbs[rowSums(mpa.org.sbs) > 0,]

OTU.df <- as.data.frame(t(mpa.org.sbs))

write.table(cbind(OTU_ID = rownames(OTU.df), OTU.df), "data/OTU.songbird", sep = "\t", quote = F, row.names = F)  

# NMDS
nmds <- metaMDS(mpa.org.sbs)
nmds.p <- as.data.frame(nmds$points)
nmds.p <- merge(metadata_ibd[c(1,2,3,4,7)], cbind(rownames(nmds.p), nmds.p), by = 1)
nmds.p$disease[nmds.p$disease == "HEALTHY"] <- "DONOR"

nmds.p$disease <- factor(nmds.p$disease, levels = c("CD", "UC", "IBS", "DONOR")) 
nmds.p$time_point[nmds.p$status == "DONOR"] <- 5

nmds_plot <- ggplot()+
    geom_point(nmds.p, mapping = aes(MDS1, MDS2, col = time_point, group = source_id, shape = disease), size = 2.5)+
    theme_bw()+
    stat_ellipse(nmds.p[nmds.p$status == "DONOR",], mapping = aes(MDS1, MDS2))+
    scale_colour_gradient(
        low = "tan1",
        high = "darkred",
        space = "Lab",
        na.value = "grey50",
        guide = "colourbar",
        aesthetics = "colour"
    )+
    scale_shape_manual(values = c(15,16,17,8))

svg(filename="figures/nmds_plot_non_sorting.svg", width=4.35, height=3.35)
nmds_plot
dev.off()

# alpha-diversity (shannon index)
div <- diversity(mpa.org.sbs)
div <- as.data.frame(div)
colnames(div) <- "shannon"
div <- merge(metadata_ibd[c(1,2,3,7)][1:65,], cbind(rownames(div), div), by = 1)
div$disease <- factor(div$disease, levels = c("CD", "UC", "IBS"))

shannon_plot <- ggplot()+
    geom_line(div, mapping = aes(time_point, shannon, group = source_id, col = disease), size = 1, alpha = 0.33)+
    stat_smooth(div, mapping = aes(time_point, shannon, col = disease), se = F)+
    facet_wrap(~disease, ncol = 3)+
    ylab("Shannon index")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

svg(filename="figures/shannon_plot.svg", width=5, height=1.75)
shannon_plot
dev.off()

# calculate Bray-Curtis dissimilarity
bray.dist <- vegdist(mpa.org.sbs)
bray.dist <- as.data.frame(as.matrix(bray.dist))

# from baseline
source_id <- unique(metadata_ibd$source_id[metadata_ibd$status != "DONOR"])
source_id <- paste0(source_id, "_")

from_baseline <- NULL
for (k in source_id){
    bray.dist.sbs <- bray.dist[which(!is.na(str_extract(colnames(bray.dist), k)))]
    bray.dist.sbs <- bray.dist.sbs[colnames(bray.dist.sbs),]
    bray.dist.sbs <- bray.dist.sbs[1]
    bray.dist.sbs <- cbind(sample_id = rownames(bray.dist.sbs), bray.dist.sbs)
    colnames(bray.dist.sbs)[2] <- "bray_dist"
    
    from_baseline <- rbind(from_baseline, bray.dist.sbs)
}

from_baseline <- merge(metadata_ibd[c(1,2,3,7)], from_baseline, by = 1)
from_baseline$disease <- factor(from_baseline$disease, levels = c("CD", "UC", "IBS"))

bray_baseline <- ggplot()+
    geom_line(from_baseline, mapping = aes(time_point, bray_dist, group = source_id, col = disease), size = 1, alpha = 0.33)+
    stat_smooth(from_baseline, mapping = aes(time_point, bray_dist, col = disease), se = F)+
    facet_wrap(~disease, ncol = 3)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

svg(filename="figures/bray_baseline.svg", width=5, height=1.75)
bray_baseline
dev.off()

# from donor
donor_id <- metadata_ibd$df_sample[metadata_ibd$status == "DONOR"]

from_donor <- NULL
for (k in donor_id){
    recip_id <- metadata_ibd$df_sample[metadata_ibd$donor_id_manual == k]
    recip_id <- recip_id[!is.na(recip_id)]
    
    bray.dist.sbs <- bray.dist[recip_id,][k]
    bray.dist.sbs <- cbind(sample_id = rownames(bray.dist.sbs), bray.dist.sbs)
    colnames(bray.dist.sbs)[2] <- "bray_dist"
    
    from_donor <- rbind(from_donor, bray.dist.sbs)
}

from_donor <- from_donor[!is.na(from_donor$bray_dist),]

from_donor <- merge(metadata_ibd[c(1,2,3,7)], from_donor, by = 1)
from_donor$disease <- factor(from_donor$disease, levels = c("CD", "UC", "IBS", "Healthy"))

bray.dist[which(!is.na(str_extract(rownames(bray.dist), "UC7")))]

bray_donor <- ggplot()+
    geom_line(from_donor, mapping = aes(time_point, bray_dist, group = source_id, col = disease), size = 1, alpha = 0.33)+
    stat_smooth(from_donor, mapping = aes(time_point, bray_dist, col = disease), se = F)+
    facet_wrap(~disease, ncol = 4)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    ylim(0,1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

svg(filename="figures/bray_donor.svg", width=5, height=1.75)
bray_donor
dev.off()