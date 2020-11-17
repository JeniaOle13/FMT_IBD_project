# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import libraries
library(ggplot2)
library(pheatmap)
library(vegan)
library(stringr)

# import metadata
metadata_ibd <- read.csv("data/metadata_ibd.txt", header = T, stringsAsFactors = F)

# import songbird results
differentials <- read.csv("output/differentials.tsv", sep = "\t")

sample_plot_data <- read.csv("output/sample_plot_data.tsv", sep = "\t")
sample_plot_data <- merge(metadata_ibd[c(1,2,3)], sample_plot_data, by = 1)
sample_plot_data$disease <- factor(sample_plot_data$disease, levels = c("CD", "UC", "IBS"))

# make quro plot
quro_plot <- ggplot()+
    geom_line(sample_plot_data, mapping = aes(time_point, Current_Natural_Log_Ratio, group = source_id, col = disease), size = 1.2, alpha = 0.33)+
    stat_smooth(sample_plot_data, mapping = aes(time_point, Current_Natural_Log_Ratio, col = disease), se = F)+
    facet_wrap(~disease)+
    ylab("Current Natural Log Ratio")+
    xlab("time point")+
    theme_bw()+
    scale_color_brewer(palette="Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

svg(filename="figures/quro_plot.svg", width=5, height=1.75)
quro_plot
dev.off()

# make differentials plot
differentials$status[differentials$effect_size > 0.75] <- "increase"
differentials$status[differentials$effect_size < -0.75] <- "dicrease"

differentials.sbs <- differentials[!is.na(differentials$status),]
differentials.sbs <- differentials.sbs[order(differentials.sbs$effect_size, decreasing = T),]

effect_size_plot <- ggplot(differentials.sbs, aes(effect_size, reorder(featureid, effect_size), fill = status))+
    geom_bar(stat="identity", width = 0.75)+
    theme_classic()+
    scale_fill_brewer(palette="Set1")+
    xlab("Effect size")+
    ylab("Species")+
    theme(legend.position = "none")+
    scale_y_discrete(position = "right")

svg(filename="figures/effect_size_plot.svg", width=4.5, height=3.5)
effect_size_plot
dev.off()