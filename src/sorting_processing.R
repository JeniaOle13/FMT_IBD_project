# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import libraries
library(stringr)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(tidyr)

# import metadat file
metadata <- read.csv("output/metadata_songbird.txt", sep = "\t")

# import and filtering mpa data
mpa.df <- read.csv("data/mpa.sorting.txt", sep = "\t", row.names = 1)
mpa.df.sbs <- mpa.df[which(!is.na(str_extract(colnames(mpa.df), "s__")))]
colnames(mpa.df.sbs) <- sapply(str_split(colnames(mpa.df.sbs), "s__"), function(x) x[2])

num.zeros <- colSums(mpa.df.sbs==0)
save.zeros  <-  (num.zeros<round(nrow(mpa.df.sbs)*0.99))
mpa.df.sbs  <-  mpa.df.sbs[,(save.zeros==TRUE)]
mpa.df.sbs <- mpa.df.sbs[rowSums(mpa.df.sbs) > 0,]

# make additional metadata file
SAMPLES <- gsub("_came_from_|_came_|both|itself|_came_from_donor|before|_settle|_stay|_not_settle|_gone", "", rownames(mpa.df.sbs))
SAMPLES <- sapply(str_split(SAMPLES, "_donor_"), function(x) x[2])
BASKET <- gsub(paste0(1:9, collapse = "_|"), "", sapply(str_split(rownames(mpa.df.sbs), "B"), function(x) tail(x, 1)))
BASKET <- gsub("came_", "", BASKET)
BASKET[BASKET == "from_before"] <- "from_baseline"

group.df <- data.frame(sample_id = rownames(mpa.df.sbs), ns_sample_id = SAMPLES, basket = BASKET, stringsAsFactors = F)
group.df.sbs <- group.df[group.df$basket %in% c("settle", "not_settle", "stay", "gone"),]

# make sorting mpa NMDS plot
mds <- metaMDS(mpa.df.sbs[rownames(mpa.df.sbs) %in% group.df.sbs$sample_id,])

mds.points <- as.data.frame(mds$points)
mds.points <- merge(group.df.sbs, cbind(sample_id = rownames(mds.points), mds.points), by = "sample_id")
mds.points$basket <- as.character(mds.points$basket)
mds.points$basket <- factor(mds.points$basket, levels = c("settle", "not_settle", "stay", "gone"))

nmds_sorting_mpa <- ggplot(mds.points, aes(MDS1, MDS2, col = basket, shape = basket))+
    geom_point(size = 2.5, stroke = 1.05)+
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    scale_shape_manual(values = c(19,1, 19,1))

svg(filename="figures/nmds_sorting_mpa.svg", width=5, height=3.75, pointsize=12)
nmds_sorting_mpa
dev.off()

# make alpha diversity plots by time
div.df <- diversity(mpa.df.sbs[rownames(mpa.df.sbs) %in% group.df.sbs$sample_id,])
div.df <- data.frame(sample_id = names(div.df), shannon = div.df)
div.df <- merge(group.df, div.df, by = 1)
div.df <- merge(metadata[c(1:3,10)], div.df[-1], by = 1)

div.df$disease <- as.character(div.df$disease)
div.df$disease <- factor(div.df$disease, levels = c("CD", "UC", "IBS"))

div.df$basket <- as.character(div.df$basket)
div.df$basket <- factor(div.df$basket, levels = c("settle", "not_settle", "stay", "gone", "from_donor", "from_both", "from_baseline", "itself"))

shannon_basket_plot <- ggplot(div.df, aes(time_point, shannon, col = basket))+
    stat_smooth(se = F, size = 1.25)+
    facet_wrap(~disease)+
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("Shannon index")+
    xlab("Time point")+
    theme(legend.position = "bottom")

# make N reads plot by time
nreads.df <- rowSums(mpa.df.sbs[rownames(mpa.df.sbs) %in% group.df.sbs$sample_id,])
nreads.df <- data.frame(sample_id = names(nreads.df), nreads = nreads.df)
nreads.df <- merge(group.df, nreads.df, by = 1)
nreads.df <- merge(metadata[c(1:3,10)], nreads.df[-1], by = 1)

nreads.df$disease <- as.character(nreads.df$disease)
nreads.df$disease <- factor(nreads.df$disease, levels = c("CD", "UC", "IBS"))

nreads.df$basket <- as.character(nreads.df$basket)
nreads.df$basket <- factor(nreads.df$basket, levels = c("settle", "not_settle", "stay", "gone", "from_donor", "from_both", "from_baseline", "itself"))

nreads_basket_plot <- ggplot(nreads.df, aes(time_point, log(nreads), col = basket))+
    stat_smooth(se = F, size = 1.25)+
    facet_wrap(~disease)+
    theme_bw()+
    scale_color_brewer(palette = "Set1")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    ylab("log N reads")+
    xlab("Time point")+
    theme(legend.position = "bottom")

#
mpa.sns <- mpa.df.sbs[rownames(mpa.df.sbs) %in% group.df.sbs$sample_id[group.df.sbs$basket %in% c("settle", "not_settle")],]
mpa.sns <- merge(group.df.sbs, cbind(settle = rownames(mpa.sns), mpa.sns), by = 1)[-1]
mpa.sns <- melt(mpa.sns)
mpa.sns <- spread(mpa.sns, basket, value, fill = 0)
mpa.sns <- mpa.sns[(mpa.sns$settle + mpa.sns$not_settle) > 0,]

mpa.sns$settle <- mpa.sns$settle+1
mpa.sns$not_settle <- mpa.sns$not_settle+1

mpa.sns$index <- mpa.sns$settle/(mpa.sns$settle+mpa.sns$not_settle)
mpa.sns$index_2[mpa.sns$index > 0.5] <- 1
mpa.sns$index_2[mpa.sns$index < 0.5] <- 0

colonization_rate_overall <- round(length(mpa.sns$index_2[mpa.sns$index_2 == 1])/(length(mpa.sns$index_2[mpa.sns$index_2 == 1])+
                                                                                      length(mpa.sns$index_2[mpa.sns$index_2 == 0])), 2)

# 
mpa.sg <- mpa.df.sbs[rownames(mpa.df.sbs) %in% group.df.sbs$sample_id[group.df.sbs$basket %in% c("stay", "gone")],]
mpa.sg <- merge(group.df.sbs, cbind(sample = rownames(mpa.sg), mpa.sg), by = 1)[-1]
mpa.sg <- melt(mpa.sg)
mpa.sg <- spread(mpa.sg, basket, value, fill = 0)
mpa.sg <- mpa.sg[mpa.sg$stay > 0 | mpa.sg$gone > 0,]

mpa.sg$stay <- mpa.sg$stay+1
mpa.sg$gone <- mpa.sg$gone+1

mpa.sg$index <- mpa.sg$stay/(mpa.sg$stay+mpa.sg$gone)
mpa.sg$index_2[mpa.sg$index > 0.5] <- 1
mpa.sg$index_2[mpa.sg$index < 0.5] <- 0

remaing_rate_overall <- round(length(mpa.sg$index_2[mpa.sg$index_2 == 1])/(length(mpa.sg$index_2[mpa.sg$index_2 == 1])+
                                                                               length(mpa.sg$index_2[mpa.sg$index_2 == 0])), 2)

density_norm_plot <- ggplot()+
    geom_density(mpa.sns, mapping = aes(index_2), col = "white", fill = 'red', alpha = 0.35)+
    geom_density(mpa.sg, mapping = aes(index_2), col = "white", fill = 'blue', alpha = 0.35)+
    theme_classic()+
    scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1.0), limits = c(-0.2, 1.2))

# rates
CD_not_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 0 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "CD"]])
CD_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 1 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "CD"]])     
UC_not_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 0 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "UC"]])
UC_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 1 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "UC"]])     
IBS_not_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 0 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "IBS"]])
IBS_settle <- length(mpa.sns$index_2[mpa.sns$index_2 == 1 & mpa.sns$ns_sample_id %in% metadata$sample_id[metadata$disease == "IBS"]])     

colonization_rate_CD <- round(CD_settle/(CD_settle+CD_not_settle), 2)
colonization_rate_UC <- round(UC_settle/(UC_settle+UC_not_settle), 2)
colonization_rate_IBS <- round(IBS_settle/(IBS_settle+IBS_not_settle), 2)

CD_gone <- length(mpa.sg$index_2[mpa.sg$index_2 == 0 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "CD"]])
CD_stay <- length(mpa.sg$index_2[mpa.sg$index_2 == 1 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "CD"]])
UC_gone <- length(mpa.sg$index_2[mpa.sg$index_2 == 0 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "UC"]])
UC_stay <- length(mpa.sg$index_2[mpa.sg$index_2 == 1 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "UC"]])
IBS_gone <- length(mpa.sg$index_2[mpa.sg$index_2 == 0 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "IBS"]])
IBS_stay <- length(mpa.sg$index_2[mpa.sg$index_2 == 1 & mpa.sg$ns_sample_id %in% metadata$sample_id[metadata$disease == "IBS"]])

remaing_rate_CD <- round(CD_stay/(CD_stay+CD_gone), 2)
remaing_rate_UC <- round(UC_stay/(UC_stay+UC_gone), 2)
remaing_rate_IBS <- round(UC_stay/(IBS_stay+IBS_gone), 2)

colonize <- c(colonization_rate_CD, colonization_rate_UC, colonization_rate_IBS)
remaing <- c(remaing_rate_CD, remaing_rate_UC, remaing_rate_IBS)
RATES <- data.frame(colonize, remaing)
rownames(RATES) <- c("CD", "UC", "IBS")

pheatmap(RATES, cluster_cols = F, display_numbers = T, fontsize = 15, color = c("white", "gray"), legend = F)
