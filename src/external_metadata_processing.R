# set work dir
workDir <- "/home/acari/github/FMT_IBD_project/"
setwd(workDir)

# Franzosa 2018
Franzosa_2018_1 <- read.csv("data/third_party_data/Franzosa_2018.csv")
colnames(Franzosa_2018_1)[2] <- "Sample"
Franzosa_2018_2 <- read.csv("data/third_party_data/Franzosa_2018_2.csv", sep = "\t", header = F)
Franzosa_2018_2 <- as.data.frame(t(Franzosa_2018_2))
rownames(Franzosa_2018_2) <- gsub("V", "", rownames(Franzosa_2018_2))
colnames(Franzosa_2018_2) <- c("Sample", "SRA_accasion", "Age", "Diagnosis")
Franzosa_2018_2 <- Franzosa_2018_2[-1,]
Franzosa_2018_2$Sample <- sub("Validation\\|", "", Franzosa_2018_2$Sample)
Franzosa_2018_2$Sample <- sub("\\|", "_", Franzosa_2018_2$Sample)

Franzosa_2018 <- merge(Franzosa_2018_1, Franzosa_2018_2, by = "Sample")[c(2,5)]
Franzosa_2018$Dataset <- "Franzosa_2019"

# He 2017
He_2017 <- read.table("data/third_party_data/He_2017.txt", header = T, row.names = 1)
colnames(He_2017) <- c("Run", "Diagnosis")
He_2017$Dataset <- "He_2017"

# 
Schirmer_2018_1 <- read.csv("data/third_party_data/Schirmer_2018.csv")[c(2,1)]
Schirmer_2018_2 <- read.csv("data/third_party_data/Schirmer_2018_2.csv")
Schirmer_2018 <- merge(Schirmer_2018_1, Schirmer_2018_2, by = 1)[c(2,9)]
Schirmer_2018$Dataset <- "Schirmer_2018"
colnames(Schirmer_2018)[2] <- "Diagnosis"

external_metadata <- rbind(Franzosa_2018, He_2017, Schirmer_2018)
external_metadata$Diagnosis <- as.character(external_metadata$Diagnosis)
external_metadata$Diagnosis[external_metadata$Diagnosis == "nonIBD"] <- "Control"

write.table(external_metadata, "data/third_party_data/external_metadata.txt", sep = "\t", quote = F, row.names = F)
