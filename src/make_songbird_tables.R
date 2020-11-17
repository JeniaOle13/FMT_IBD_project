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

# read metadata
metadata_songbird <- read.table("output/metadata_songbird.txt", header = T, stringsAsFactors = F, sep = "\t")

# make extended metadata file
DONOR <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[1])
BASKET  <- sapply(str_split(rownames(mpa.org.sbs), "_donor_"), function(x) x[2])
SAMPLE <- gsub("_came_from_before|_came_from_donor|_came_from_both|_came_itself|_settle|_not_settle|_stay|_gone", "", BASKET)
BASKET <- sapply(str_split(BASKET, "T"), function(x) x[2])
BASKET <- gsub(paste0(paste0(1:9, "_"), collapse = "|"), "", sapply(str_split(BASKET, "B"), function(x) x[2]))

basket_meta <- merge(metadata_ibd[c(1:3,5,7)], data.frame(SAMPLE, full_name = rownames(mpa.org.sbs), basket = BASKET), by = 1)
basket_meta <- basket_meta[c(6,1:5,7)]

# settle/not_settle
df.settle <- mpa.org.sbs[which(!is.na(str_extract(rownames(mpa.org.sbs), "settle"))),]
df.settle <- df.settle[colSums(df.settle) > 0]
df.settle <- df.settle[rowSums(df.settle) > 0,]
df.settle <- df.settle[order(rownames(df.settle)),]

meta_settle <- basket_meta[as.character(basket_meta$full_name) %in% rownames(df.settle),]
colnames(metadata_songbird)[1] <- "df_sample"
meta_settle <- merge(metadata_songbird[c(1,4:5,7:9)], meta_settle, by = "df_sample")
meta_settle <- meta_settle[c(7,1:6,8:12)]
meta_settle$full_name <- as.character(meta_settle$full_name)
meta_settle <- meta_settle[order(meta_settle$full_name),]
colnames(meta_settle)[1] <- "sample_id"

# write.table(t(df.settle), "output/OTU.songbird.settle", quote = F, sep = "\t")
# write.table(meta_settle, "output/metadata_songbird_settle.txt", quote = F, sep = "\t", row.names = F)

# stay/gone
df.stay <- mpa.org.sbs[which(!is.na(str_extract(rownames(mpa.org.sbs), "stay|gone"))),]
df.stay <- df.stay[colSums(df.stay) > 0]
df.stay <- df.stay[rowSums(df.stay) > 0,]
df.stay <- df.stay[order(rownames(df.stay)),]

meta_stay <- basket_meta[as.character(basket_meta$full_name) %in% rownames(df.stay),]
colnames(metadata_songbird)[1] <- "df_sample"
meta_stay <- merge(metadata_songbird[c(1,4:5,7:9)], meta_stay, by = "df_sample")
meta_stay <- meta_stay[c(7,1:6,8:12)]
meta_stay$full_name <- as.character(meta_stay$full_name)
meta_stay <- meta_stay[order(meta_stay$full_name),]
colnames(meta_stay)[1] <- "sample_id"

# write.table(t(df.stay), "output/OTU.songbird.stay", quote = F, sep = "\t")
# write.table(meta_stay, "output/metadata_songbird_stay.txt", quote = F, sep = "\t", row.names = F)

df.donor <- mpa.org.sbs[which(!is.na(str_extract(rownames(mpa.org.sbs), "settle"))),]
df.donor <- df.donor[colSums(df.donor) > 0]
df.donor <- df.donor[rowSums(df.donor) > 0,]
df.donor <- df.donor[order(rownames(df.donor)),]

colSums(mpa.sorting.org==0)