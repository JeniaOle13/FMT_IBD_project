# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# import libraries
library(stringr)

# import data
donor_info <- read.csv("data/metadata/donor_info.tsv", sep = "\t", stringsAsFactors = F)[c(1,4)]

clin_index <- read.csv("data/metadata/clin_index.tsv", sep = "\t")

df_samples_meta_final <- read.csv("data/metadata/df_samples_meta_final.tsv", sep = "\t", stringsAsFactors = F)[c(1,2,4)]
df_samples_meta_final <- df_samples_meta_final[df_samples_meta_final$representative == "ok",]
df_samples_meta_final <- df_samples_meta_final[-3]

link_donor <- merge(df_samples_meta_final, donor_info, by = "source_id")
link_donor$time_point <- as.numeric(sub("T", "", sapply(str_split(link_donor$df_sample, "_"), function(x) x[2])))

link_donor <- merge(link_donor, clin_index, by = c("source_id", "time_point"), all = T)[-6]
link_donor <- link_donor[!is.na(link_donor$df_sample),]
link_donor <- link_donor[c(3,1,4,5,2)]

write.table(link_donor, "data/metadata.txt", quote = F, row.names = F, sep = "\t")

# import data
df <- read.csv("data/metadata/df.non_sorting.ibd.mpa.out", sep = "\t", row.names = 1)

df.sbs <- df[which(!is.na(str_extract(colnames(df), "s__")))]
colnames(df.sbs) <- gsub("s__", "", sapply(str_split(colnames(df.sbs), "\\."), function(x) x[7]))
df.sbs <- df.sbs[order(colSums(df.sbs), decreasing = T)]

# filtration 
Numzz <- colSums(df.sbs==0)
Savezz <-  (Numzz<nrow(df.sbs)*0.75)
df.sbs <- df.sbs[Savezz]

df_samples_meta_final <- read.csv("data/metadata/df_samples_meta_final.tsv", sep = "\t", stringsAsFactors = F)[c(1,2,4)]
df_samples_meta_final <- df_samples_meta_final[df_samples_meta_final$representative == "ok",]
df.sbs <- df.sbs[rownames(df.sbs) %in% c(df_samples_meta_final$df_sample, "SFM028_T1_B1"),]

donor_ids <- rownames(df.sbs)[!rownames(df.sbs) %in% metadata$df_sample] 
donor_ids_delete <- gsub(paste0(paste0("_B", c(1:9)), collapse = "|"), "", donor_ids)
donor_ids_delete <- donor_ids_delete[!donor_ids_delete %in% unique(metadata$donor_id_manual)]

df.sbs <- df.sbs[!rownames(df.sbs) %in% donor_ids[c(1,10)],]

write.table(df.sbs, "data/mpa_ibd_non_sorting.txt", row.names = F, quote = F, sep = "\t")

## 

setwd("/home/acari/github/FMT_IBD_project/data/metadata/strain_finder")

samples.meta_2 <- read.csv("samples.meta_2.csv", sep = "\t", stringsAsFactors = F) 
sra_link <- read.csv("sra_link.csv")[c(2,1)]
metadata_cdi <- read.csv("metadata_cdi.csv", stringsAsFactors = F)

LIST <- str_split(samples.meta_2$gids, "\\,")
names(LIST) <- samples.meta_2$sample

LIST_1 <- sapply(LIST, function(x) x[1])
LIST_2 <- sapply(LIST, function(x) x[2])
LIST_2 <- LIST_2[!is.na(LIST_2)]

LIST <- c(LIST_1, LIST_2)

df <- data.frame(sample_id = names(LIST), gids = LIST, stringsAsFactors = F)
df <- df[order(df$sample_id, decreasing = F),]
df <- df[c(2,1)]

df <- merge(sra_link, df, by = 1)

metadata_cdi <- merge(metadata_cdi, df[c(3,1,2)], by = 1, all = T)

write.table(metadata_cdi, "/home/acari/github/FMT_IBD_project/data/metadata_cdi.txt", sep = "\t", quote = F)

metadata_cdi <- read.csv("/home/acari/github/FMT_IBD_project/data/metadata_cdi.txt", stringsAsFactors = F)

write.csv(df[order(df$sample_id),], "df.all")

# 
# import data
df <- read.csv("data/metadata/strain_finder/df.non_sorting.cdi.mpa.out", sep = "\t", row.names = 1)
additional_meta <- read.csv("data/additional_meta_cdi.txt", stringsAsFactors = F)

df.sbs <- df[which(!is.na(str_extract(colnames(df), "s__")))]
colnames(df.sbs) <- gsub("s__", "", sapply(str_split(colnames(df.sbs), "\\."), function(x) x[7]))
df.sbs <- df.sbs[order(colSums(df.sbs), decreasing = T)]

# filtration 
Numzz <- colSums(df.sbs==0)
Savezz <-  (Numzz<nrow(df.sbs)*0.75)
df.sbs <- df.sbs[Savezz]

df.sbs <- df.sbs[rownames(df.sbs) %in% additional_meta$run,]

write.table(df.sbs, "data/mpa_cdi_non_sorting.txt", row.names = F, quote = F, sep = "\t")
