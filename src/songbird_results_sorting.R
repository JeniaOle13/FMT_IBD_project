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

# import species data
mpa.sorting <- read.table("data/mpa.sorting.txt", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)
mpa.sorting.org <- mpa.sorting[which(!is.na(str_extract(colnames(mpa.sorting), "s__")))]
colnames(mpa.sorting.org) <- sub("s__", "", sapply(str_split(colnames(mpa.sorting.org), "\\."), function(x) x[7]))

# filtering data
Numzz <- colSums(mpa.sorting.org==0)
Savezz <-  (Numzz<nrow(mpa.sorting.org)*0.95)

mpa.org.sbs <- mpa.sorting.org[Savezz]
mpa.org.sbs <- mpa.org.sbs[order(colSums(mpa.org.sbs), decreasing = T)]

ids <- c(metadata_ibd$df_sample, metadata_ibd$donor_id_manual)
ids <- unique(ids[!is.na(ids)])

mpa.org.sbs <- mpa.org.sbs[colSums(mpa.org.sbs) > 10000]

# import songbird results
differentials_settle <- read.csv("output/differentials_settle.tsv", sep = "\t")[c(1,38)]
differentials_stay <- read.csv("output/differentials_stay.tsv", sep = "\t")[c(1,38)]

colnames(differentials_settle) <- c("species", "settle")
colnames(differentials_stay) <- c("species", "stay")

differentials_settle <- differentials_settle[differentials_settle$settle > 3 | differentials_settle$settle < -3,]
differentials_stay <- differentials_stay[differentials_stay$stay > 3 | differentials_stay$stay < -3,]

differentials <- merge(differentials_settle, differentials_stay, by = 1)

# settle_up <- differentials_settle[differentials_settle$settle > 3,]
# settle_up <- settle_up[order(settle_up$settle, decreasing = T),][c(1:10),]
# 
# settle_down <- differentials_settle[differentials_settle$settle < -3,]
# settle_down <- settle_down[order(settle_down$settle),][c(1:10),]
# 
# differentials_settle <- rbind(settle_up, settle_down)
# 
# stay_up <- differentials_stay[differentials_stay$stay > 3,]
# stay_up <- stay_up[order(stay_up$stay, decreasing = T),][c(1:10),]
# 
# stay_down <- differentials_stay[differentials_stay$stay < -3,]
# stay_down <- stay_down[order(stay_down$stay),][c(1:10),]
# 
# differentials_stay <- rbind(stay_up, stay_down)
# 
# differentials <- merge(differentials_settle, differentials_stay, by = 1, all = T)

# differentials[is.na(differentials)] <- 0

plot(differentials$settle, differentials$stay)

df.melt <- melt(differentials)
colnames(df.melt) <- c("Species", "Group", "Effect_size")
df.melt$Color[df.melt$Effect_size > 0] <- "increase"
df.melt$Color[df.melt$Effect_size < 0] <- "decrease"


ggplot(df.melt, aes(Effect_size, reorder(Species, Effect_size), fill = Color))+
    geom_bar(stat = "identity", width = 0.75)+
    facet_wrap(~Group)+
    theme_classic()+
    xlab("Species")+
    ylab("Effect size")+
    scale_fill_brewer(palette="Set1")+
    theme(legend.position = "none")+
    scale_y_discrete(position = "right")

# make subset biplot
# subset_species <- unique(c(as.character(differentials_settle$species), as.character(differentials_stay$species)))
subset_species <- unique(as.character(df.melt$Species))

mpa.org.sbs <- mpa.org.sbs[subset_species]

# make extended metadata file
DONOR <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[1])
BASKET  <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[2])
SAMPLE <- gsub("_came_from_before|_came_from_donor|_came_from_both|_came_itself|_settle|_not_settle|_stay|_gone", "", BASKET)
BASKET <- sapply(str_split(BASKET, "T"), function(x) x[2])
BASKET <- gsub(paste0(paste0(1:9, "_"), collapse = "|"), "", sapply(str_split(BASKET, "B"), function(x) x[2]))

basket_meta <- merge(metadata_ibd[c(1:3,5,7)], data.frame(SAMPLE, full_name = rownames(mpa.org.sbs), basket = BASKET), by = 1)
basket_meta <- basket_meta[c(6,1:5,7)]
basket_meta$n_reads <- rowSums(mpa.org.sbs)

basket_meta.sbs <- basket_meta[which(is.na(str_extract(basket_meta$basket, "came") )),]
basket_meta.sbs$basket <- as.character(basket_meta.sbs$basket)
basket_meta.sbs$full_name <- as.character(basket_meta.sbs$full_name)

# make biplot
mpa.org.sbs_2 <- mpa.org.sbs[rownames(mpa.org.sbs) %in% basket_meta.sbs$full_name,]
mpa.org.sbs_2 <- mpa.org.sbs_2[rowSums(mpa.org.sbs_2) > 0,]

mds <- metaMDS(mpa.org.sbs_2)
mds.points <- as.data.frame(mds$points)

mds.points <- merge(basket_meta.sbs, cbind(rownames(mds.points), mds.points), by = 1)
mds.points$basket <- factor(mds.points$basket, levels = c("settle", "not_settle", "stay", "gone"))
mds.points$disease <- factor(mds.points$disease, levels = c("CD", "UC", "IBS"))

nmds_sorting_subset_plot <- ggplot(mds.points, aes(MDS1, MDS2, shape = disease, col = basket))+
    geom_point(size = 2.5, alpha = 0.85)+
    scale_color_brewer(palette = "Set1")+
    theme_bw()+
    scale_shape_manual(values = c(15,16,17))

svg(filename="figures/nmds_sorting_subset_plot.svg", width=4.35, height=3.35)
nmds_sorting_subset_plot
dev.off()