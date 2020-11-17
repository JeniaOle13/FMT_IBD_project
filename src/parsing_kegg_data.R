setwd("/home/acari/github/FMT_IBD_project/")

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

DEFINITION <- read.table("data/kegg/DEFINITION")
DEFINITION$V1 <- sub("DEFINITION  ", "", DEFINITION$V1)

ENTRY <- read.table("data/kegg/ENTRY", sep = "\t")
ENTRY$V1 <- sapply(str_split(ENTRY$V1, " "), function(x) x[8])

ko.df <- cbind(ENTRY, DEFINITION)
colnames(ko.df) <- c("ko_id", "ko_name")

write.table(ko.df, "data/kegg/ko.df.txt", quote = F, sep = "\t", row.names = F)

ko_pathway <- read.table("data/kegg/ko_pathway.list", sep = "\t", header = F, stringsAsFactors = F)
ko_pathway$V1 <- sub("ko:", "", ko_pathway$V1)
ko_pathway$V2 <- sub("path:", "", ko_pathway$V2)
ko_pathway <- ko_pathway[which(is.na(str_extract(ko_pathway$V2, "map"))),]

write.table(ko_pathway, "data/kegg/ko_pathway.list", quote = F, sep = "\t", row.names = F)
