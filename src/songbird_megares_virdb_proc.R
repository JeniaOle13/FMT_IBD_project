# import libraries
library(dplyr)
library(stringr)
library(tidyr)
library(pheatmap)

# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import data
megares <- read.csv("data/megares.out.txt", sep = "\t", stringsAsFactors = F)
metadata_ibd <- read.csv("data/metadata_ibd.txt", stringsAsFactors = F)
metadata_songbird <- read.csv("output/metadata_songbird.txt", stringsAsFactors = F, sep = "\t")

virdb_setA_depth <- read.csv("data/virdb_setA_depth.txt", sep = "\t", stringsAsFactors = F)
virdb_setA_idxstats <- read.csv("data/virdb_setA_idxstats.txt", sep = "\t", stringsAsFactors = F)
virdb_setB_depth <- read.csv("data/virdb_setB_depth.txt", sep = "\t", stringsAsFactors = F)
virdb_setB_idxstats <- read.csv("data/virdb_setB_idxstats.txt", sep = "\t", stringsAsFactors = F)

# make megares df
gene_id <- sapply(str_split(megares$Gene, "\\|"), function(x) x[1])
superfamily_id <- sapply(str_split(megares$Gene, "\\|"), function(x) x[2])
family_id <- sapply(str_split(megares$Gene, "\\|"), function(x) x[3])
mechanism_id <- sapply(str_split(megares$Gene, "\\|"), function(x) x[4])
group_id <- sapply(str_split(megares$Gene, "\\|"), function(x) x[5])

gene_info <- data.frame(gene_id, superfamily_id, family_id, mechanism_id, group_id, stringsAsFactors = F)
gene_info <- unique(gene_info)

megares$Gene <- group_id
megares <- as.data.frame(megares %>% group_by(Sample, Gene) %>% summarise(Hits = sum(Hits)))

df_megares <- spread(megares, Gene, Hits, fill = 0)
rownames(df_megares) <- df_megares$Sample
df_megares <- df_megares[-1]

df_megares <- df_megares[rownames(df_megares) %in% metadata_songbird$sample_id,]
metadata_songbird_megares <- metadata_songbird[metadata_songbird$sample_id %in% rownames(df_megares),]

# filtering data
Numzz <- colSums(df_megares==0)
Savezz <-  (Numzz<nrow(df_megares)*0.75)
df_megares_sbs <- df_megares[Savezz]
df_megares_sbs <- df_megares_sbs[order(colSums(df_megares_sbs), decreasing = T)]

write.table(t(df_megares_sbs), "output/OTU.songbird.megares", quote = F, sep = "\t")
write.table(metadata_songbird_megares, "output/metadata_songbird_megares.txt", quote = F, sep = "\t", row.names = F)

# make virdb df
virdb_setB <- merge(virdb_setB_depth, virdb_setB_idxstats, by = c("gene_name", "sample"))[-3]

df_setB <- spread(virdb_setB, gene_name, hits, fill = 0)
rownames(df_setB) <- df_setB$sample
df_setB <- df_setB[-1]

df_setB <- df_setB[rownames(df_setB) %in% metadata_songbird$sample_id,]
metadata_songbird_virdb <- metadata_songbird[metadata_songbird$sample_id %in% rownames(df_setB),]

Numzz <- colSums(df_setB==0)
Savezz <-  (Numzz<nrow(df_setB)*0.75)
df_setB_sbs <- df_setB[Savezz]
df_setB_sbs <- df_setB_sbs[order(colSums(df_setB_sbs), decreasing = T)]

write.table(t(df_setB_sbs), "output/OTU.songbird.virdb", quote = F, sep = "\t")
write.table(metadata_songbird_virdb, "output/metadata_songbird_virdb.txt", quote = F, sep = "\t", row.names = F)