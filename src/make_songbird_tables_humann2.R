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
metadata_songbird <- read.table("output/metadata_songbird.txt", header = T, stringsAsFactors = F, sep = "\t")

# import data
humann2_sorting <- as.data.frame(t(read.table("data/humann2_sorting.tsv", header = T, stringsAsFactors = F, sep = "\t", row.names = 1)[-1,]))
rownames(humann2_sorting) <- sub("_Abundance.RPKs", "", rownames(humann2_sorting))

# filtering data
Numzz <- colSums(humann2_sorting==0)
Savezz <-  (Numzz<nrow(humann2_sorting)*0.75)
humann2_sorting.sbs <- humann2_sorting[Savezz]

ids <- c(metadata_ibd$df_sample, metadata_ibd$donor_id_manual)
ids <- unique(ids[!is.na(ids)])

# make extended metadata file
DONOR <- sapply(str_split(rownames(humann2_sorting.sbs), "_donor_"), function(x) x[1])
BASKET  <- sapply(str_split(rownames(humann2_sorting.sbs), "_donor_"), function(x) x[2])
SAMPLE <- gsub("_came_from_before|_came_from_donor|_came_from_both|_came_itself|_settle|_not_settle|_stay|_gone", "", BASKET)
BASKET <- sapply(str_split(BASKET, "T"), function(x) x[2])
BASKET <- gsub(paste0(paste0(1:9, "_"), collapse = "|"), "", sapply(str_split(BASKET, "B"), function(x) x[2]))

basket_meta <- merge(metadata_ibd[c(1:3,5,7)], data.frame(SAMPLE, full_name = rownames(humann2_sorting.sbs), basket = BASKET), by = 1)
basket_meta <- basket_meta[c(6,1:5,7)]

# settle/not_settle
df.settle <- humann2_sorting.sbs[which(!is.na(str_extract(rownames(humann2_sorting.sbs), "settle"))),]
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

write.table(t(df.settle), "output/humann2_settle_songbird", quote = F, sep = "\t")
write.table(meta_settle, "output/metadata_songbird_humann2_settle", quote = F, sep = "\t", row.names = F)

# stay/gone
df.stay <- humann2_sorting.sbs[which(!is.na(str_extract(rownames(humann2_sorting.sbs), "stay|gone"))),]
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

write.table(t(df.stay), "output/humann2_stay_songbird", quote = F, sep = "\t")
write.table(meta_stay, "output/metadata_songbird_humann2_stay", quote = F, sep = "\t", row.names = F)

# from_donor/from_baseline
df.donor <- humann2_sorting.sbs[which(!is.na(str_extract(rownames(humann2_sorting.sbs), "came_from_donor"))),]