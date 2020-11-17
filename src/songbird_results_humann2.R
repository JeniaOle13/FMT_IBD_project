# import libraries
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)
library(vegan)
library(ggplot2)

# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import kegg data
brite_table <- read.csv("data/kegg/brite_table.csv", stringsAsFactors = F, sep = "\t")
brite_table <- brite_table[c(3,4,1,2)]
brite_table$Pathway_id[nchar(brite_table$Pathway_id) == 3] <- paste0("ko00", brite_table$Pathway_id[nchar(brite_table$Pathway_id) == 3])
brite_table$Pathway_id[nchar(brite_table$Pathway_id) == 4] <- paste0("ko0", brite_table$Pathway_id[nchar(brite_table$Pathway_id) == 4])

DEFINITION <- read.table("data/kegg/DEFINITION")
DEFINITION$V1 <- sub("DEFINITION  ", "", DEFINITION$V1)
ENTRY <- read.table("data/kegg/ENTRY", sep = "\t")
ENTRY$V1 <- sapply(str_split(ENTRY$V1, " "), function(x) x[8])
ko_table <- cbind(ENTRY, DEFINITION)
colnames(ko_table) <- c("ko_id", "ko_name")

ko_pathway <- read.csv("data/kegg/ko_pathway.list", sep = "\t")
colnames(ko_pathway) <- c("ko", "pathway")

brite_ko <- read.csv("data/kegg/brite_ko.list", sep = "\t", header = F, stringsAsFactors = F)
colnames(brite_ko) <- c("brite", "ko")

brite_ko$brite <- sub("br:", "", brite_ko$brite)
brite_ko$ko <- sub("ko:", "", brite_ko$ko)

brite_name <- read.csv("data/kegg/brite_name_table.csv", sep = "\t", header = T, stringsAsFactors = F)

# import data
metadata_ibd <- read.csv("data/metadata_ibd.txt", stringsAsFactors = F)
metadata_songbird <- read.csv("output/metadata_songbird.txt", stringsAsFactors = F, sep = "\t")

df.kegg <- as.data.frame(t(read.csv("data/humann2_non_sorting.tsv", sep = "\t", row.names = 1)[-1,]))
rownames(df.kegg) <- gsub("_Abundance.RPKs", "", rownames(df.kegg))
df.kegg.sbs <- df.kegg[rownames(df.kegg) %in% metadata_ibd$df_sample,]

Numzz <- colSums(df.kegg.sbs==0)
Savezz <-  (Numzz<nrow(df.kegg.sbs)*0.75)
df.kegg.sbs <- df.kegg.sbs[Savezz]
df.kegg.sbs <- df.kegg.sbs[order(colSums(df.kegg.sbs), decreasing = T)]

# import songbird results
differentials_humann2 <- read.csv("output/differentials_humann2.tsv", sep = "\t" , stringsAsFactors = F)[c(1,21)]
colnames(differentials_humann2)[2] <- "effect_size"
diff_sbs <- differentials_humann2[differentials_humann2$effect_size > 1 | differentials_humann2$effect_size < -1,]

sample_plot_data <- read.csv("output/sample_plot_data_humann2.tsv", sep = "\t" , stringsAsFactors = F)[-4]
sample_plot_data <- merge(metadata_ibd[c(1,2,3)], sample_plot_data, by = 1)
sample_plot_data$disease <- factor(sample_plot_data$disease, levels = c("CD", "UC", "IBS"))

# make quro plot
quro_plot_humann2 <- ggplot()+
    geom_line(sample_plot_data, mapping = aes(time_point, Current_Natural_Log_Ratio, group = source_id, col = disease), size = 1.2, alpha = 0.33)+
    stat_smooth(sample_plot_data, mapping = aes(time_point, Current_Natural_Log_Ratio, col = disease), se = F)+
    facet_wrap(~disease)+
    ylab("Current Natural Log Ratio")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# svg(filename="figures/quro_plot_humann2.svg", width=5, height=1.75)
# quro_plot_humann2
# dev.off()

# make biplot 
# NMDS
df.kegg.sbs_2 <- df.kegg.sbs[colnames(df.kegg.sbs) %in% diff_sbs$featureid]
df.kegg.sbs_2 <- df.kegg.sbs_2[rowSums(df.kegg.sbs_2) > 0,]

nmds <- metaMDS(df.kegg.sbs_2)
nmds.p <- as.data.frame(nmds$points)
nmds.p <- merge(metadata_ibd[c(1,2,3,4,7)], cbind(rownames(nmds.p), nmds.p), by = 1)
nmds.p$disease[nmds.p$disease == "HEALTHY"] <- "DONOR"

nmds.p$disease <- factor(nmds.p$disease, levels = c("CD", "UC", "IBS", "DONOR")) 
nmds.p$time_point[nmds.p$status == "DONOR"] <- 5

nmds_plot_humann2 <- ggplot()+
    geom_point(nmds.p, mapping = aes(MDS1, MDS2, col = time_point, group = source_id, shape = disease), size = 2.5)+
    theme_bw()+
    xlim(c(-1,1.5))+
    ylim(c(-1,1.5))+
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

# svg(filename="figures/nmds_plot_humann2_non_sorting.svg", width=4.35, height=3.35)
# nmds_plot_humann2
# dev.off()

# calculate Bray-Curtis dissimilarity
bray.dist <- vegdist(df.kegg.sbs_2)
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

bray_baseline_humann2 <- ggplot()+
    geom_line(from_baseline, mapping = aes(time_point, bray_dist, group = source_id, col = disease), size = 1, alpha = 0.33)+
    stat_smooth(from_baseline, mapping = aes(time_point, bray_dist, col = disease), se = F)+
    facet_wrap(~disease, ncol = 3)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# svg(filename="figures/bray_baseline_humann2.svg", width=5, height=1.75)
# bray_baseline_humann2
# dev.off()

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

bray_donor_humann2 <- ggplot()+
    geom_line(from_donor, mapping = aes(time_point, bray_dist, group = source_id, col = disease), size = 1, alpha = 0.33)+
    stat_smooth(from_donor, mapping = aes(time_point, bray_dist, col = disease), se = F)+
    facet_wrap(~disease, ncol = 4)+
    ylab("Bray-Curtis dissimilarity")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    ylim(0,1)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# svg(filename="figures/bray_donor_humann2.svg", width=5, height=1.75)
# bray_donor_humann2
# dev.off()

# make diff density plot
diff_sbs <- diff_sbs[order(diff_sbs$effect_size, decreasing = T),]

density_plot_humann2 <- ggplot(diff_sbs, aes(effect_size))+
    geom_density(alpha = 0.75, fill = "black")+
    theme_bw()+
    theme(legend.position = "none")

# svg(filename="figures/density_plot_humann2.svg", width=5, height=1.75)
# density_plot_humann2
# dev.off()

#########
diff_sbs_2 <- merge(brite_table, merge(ko_pathway, diff_sbs, by = 1)[c(2,1,3)], by = 1)
diff_sbs_2 <- diff_sbs_2[as.character(diff_sbs_2$Pathway_id) != 'ko01100',]

UP <- diff_sbs_2$Brite_2[diff_sbs_2$effect_size > 0]
DOWN <- diff_sbs_2$Brite_2[diff_sbs_2$effect_size < 0]

UP <- as.data.frame(table(UP))
DOWN <- as.data.frame(table(DOWN))
DOWN$Freq <- DOWN$Freq*-1
colnames(DOWN)[1] <- "pathway"
colnames(UP)[1] <- "pathway"

PATHS <- rbind(DOWN, UP)
PATHS$col[PATHS$Freq > 0] <- "increase"
PATHS$col[PATHS$Freq < 0] <- "decrease"
PATHS$pathway <- as.character(PATHS$pathway)

PATHS.sbs <- PATHS[abs(PATHS$Freq) > 1,]
 
brite_humann2_plot <- ggplot(PATHS, aes(Freq, reorder(pathway, abs(Freq)), fill = col))+
    geom_bar(stat = "identity", width = 0.75)+
    # coord_flip()+
    theme_classic()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "bottom")+
    scale_fill_brewer(palette="Set1")+
    ylab("KEGG BRITE")+
    xlim(c(-45,45))+
    scale_y_discrete(position = "right")

svg(filename="figures/brite_humann2_plot.svg", width=6.0, height=3.75)
brite_humann2_plot
dev.off()

length(diff_sbs$featureid[diff_sbs$effect_size > 0])
length(diff_sbs$featureid[diff_sbs$effect_size < 0])

#########
diff_sbs_3 <- merge(brite_ko[c(2,1)], diff_sbs, by = 1)
diff_sbs_3 <- diff_sbs_3[!as.character(diff_sbs_3$brite) %in% c("ko00001", "ko00000", "ko00002"),]

diff_sbs_3$brite[diff_sbs_3$effect_size > 0]

UP <- diff_sbs_3$brite[diff_sbs_3$effect_size > 0]
DOWN <- diff_sbs_3$brite[diff_sbs_3$effect_size < 0]

UP <- as.data.frame(table(UP))
DOWN <- as.data.frame(table(DOWN))
DOWN$Freq <- DOWN$Freq*-1

colnames(DOWN)[1] <- "Brite"
colnames(UP)[1] <- "Brite"

PATHS <- rbind(DOWN, UP)
PATHS$col[PATHS$Freq > 0] <- "increase"
PATHS$col[PATHS$Freq < 0] <- "decrease"
PATHS$Brite <- as.character(PATHS$Brite)

PATHS.sbs <- PATHS[abs(PATHS$Freq) > 1,]
PATHS.sbs <- merge(brite_name, PATHS.sbs, by = 1)
PATHS.sbs$Brite <- paste(PATHS.sbs$brite_id, PATHS.sbs$brite_name, sep = "| ")

brite_humann2_plot <- ggplot(PATHS.sbs, aes(Freq, reorder(Brite, abs(Freq)), fill = col))+
    geom_bar(stat = "identity", width = 0.75)+
    # coord_flip()+
    theme_classic()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "bottom")+
    scale_fill_brewer(palette="Set1")+
    ylab("KEGG BRITE")+
    xlim(c(-75,75))+
    scale_y_discrete(position = "right")

svg(filename="figures/brite_humann2_plot.svg", width=7.5, height=3.75)
brite_humann2_plot
dev.off()

ko_name <- merge(ko_table, diff_sbs, by = 1)
ko_name <- ko_name[order(ko_name$effect_size, decreasing = T),]

ko_name <- merge(ko_name, ko_pathway, by = 1)

UP <- ko_name$pathway[ko_name$effect_size > 0]
DOWN <- ko_name$pathway[ko_name$effect_size < 0]

UP <- as.data.frame(table(as.character(UP)))
DOWN <- as.data.frame(table(as.character(DOWN)))
DOWN$Freq <- DOWN$Freq*-1

colnames(DOWN)[1] <- "Brite"
colnames(UP)[1] <- "Brite"

PATHS <- rbind(DOWN, UP)
PATHS$col[PATHS$Freq > 0] <- "increase"
PATHS$col[PATHS$Freq < 0] <- "decrease"
PATHS$Brite <- as.character(PATHS$Brite)

PATHS.sbs <- PATHS[abs(PATHS$Freq) > 1,]
PATHS.sbs <- merge(brite_table, PATHS.sbs, by = 1)
PATHS.sbs$Pathway <- paste(PATHS.sbs$Pathway_id, PATHS.sbs$Pathway_name, sep = "| ")

ggplot(PATHS.sbs, aes(Freq, reorder(Pathway, abs(Freq)), fill = col))+
    geom_bar(stat = "identity", width = 0.75)+
    # coord_flip()+
    theme_classic()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "bottom")+
    scale_fill_brewer(palette="Set1")+
    ylab("KEGG BRITE")+
    # xlim(c(-75,75))+
    scale_y_discrete(position = "right")