# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import libraries
library(ggplot2)
library(pheatmap)
library(vegan)
library(stringr)
library(dplyr)
library(reshape2)

# import metadata
metadata_ibd <- read.csv("data/metadata_ibd.txt", header = T, stringsAsFactors = F)

# import data
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

# make extended metadata file
DONOR <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[1])
BASKET  <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[2])
SAMPLE <- gsub("_came_from_before|_came_from_donor|_came_from_both|_came_itself|_settle|_not_settle|_stay|_gone", "", BASKET)
BASKET <- sapply(str_split(BASKET, "T"), function(x) x[2])
BASKET <- gsub(paste0(paste0(1:9, "_"), collapse = "|"), "", sapply(str_split(BASKET, "B"), function(x) x[2]))

basket_meta <- merge(metadata_ibd[c(1:3,5,7)], data.frame(SAMPLE, full_name = rownames(mpa.org.sbs), basket = BASKET), by = 1)
basket_meta <- basket_meta[c(6,1:5,7)]
basket_meta$n_reads <- rowSums(mpa.org.sbs)

basket_meta.sbs <- basket_meta[which(!is.na(str_extract(basket_meta$basket, "came") )),]

basket_meta.sbs <- merge(basket_meta.sbs, as.data.frame(basket_meta.sbs %>% group_by(df_sample) %>% 
                        summarise(sum_reads= sum(n_reads))), by = "df_sample")

basket_meta.sbs$relab <- (basket_meta.sbs$n_reads/basket_meta.sbs$sum_reads)*100
basket_meta.sbs$basket <- factor(as.character(basket_meta.sbs$basket), levels = c("came_from_donor", "came_from_both", 
                                                           "came_from_before", "came_itself"))

# barplots
basket_barplots <- ggplot(basket_meta.sbs, aes(time_point, relab, group = basket, fill = basket))+
    geom_bar(stat = "identity", width = 0.55)+
    facet_wrap(~source_id, ncol = 1)+
    coord_flip()+
    scale_fill_manual(values = c("#E41A1C", "#FF7F00", "#4DAF4A", "gray40"))+
    theme_linedraw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle=90, vjust=0.5, size=6.5))+
    theme(legend.position = "bottom")

svg(filename="figures/basket_barplots.svg", width=5, height=30)
basket_barplots
dev.off()