## This the first script of the project, dealing mainly with data preprocessing
#Note that all the data extracted or generated are stored in data folder

##Required libraries
library(dplyr)
library(GEOquery)
library(stringr)

##Reading the reads count matrix
count_raw <- read.table(file= "/Users/pragata/Codes/R_codes/Gene_exp_analysis/data/GSE81089_readcounts_featurecounts.tsv", header = TRUE, sep="\t", row.names = 1) 
head(count_raw)

##Loading the experiment metadata
gse <- getGEO("GSE81089", GSEMatrix = TRUE)
metadata<- pData(gse[[1]])
head(metadata)

##Sub-setting metadata
metadata_new <- select(metadata, c(1,8,10,11,12,14,15,17,18))

##Pre-processing metadata
#Renaming columns
metadata_new <-  rename(metadata_new, samples= title, source_tissue= source_name_ch1, condition= characteristics_ch1, stage= characteristics_ch1.1, histology= characteristics_ch1.2, age= characteristics_ch1.4, gender= characteristics_ch1.5, died= characteristics_ch1.7, Smoking=characteristics_ch1.8)

#Mutating rows
metadata_new <- metadata_new %>% 
  mutate(stage=str_replace_all(stage,"stage tnm:","")) %>% 
  mutate(gender=str_replace_all(gender,"gender:","")) %>% 
  mutate(age=str_replace_all(age,"age:","")) %>% 
  mutate(histology=str_replace_all(histology,"histology:","")) %>% 
  mutate(Smoking=str_replace_all(Smoking,"smoking:","")) %>% 
  mutate(died=str_replace_all(died,"dead:","")) %>% 
  mutate(source_tissue=str_replace_all(source_tissue,"Human","")) %>% 
  mutate(condition=str_replace_all(condition," tumor (t) or normal (n):",""))

#Replacing the meaning of indexes across various columns according to the overall design of the experiment
metadata_new <- metadata_new %>%
  mutate(
    died = str_trim(died),  
    died = case_when(
      died == "0" ~ "no",  
      died == "1" ~ "yes", 
      TRUE ~ died          
    )
  ) %>%
  mutate(
    stage = str_trim(stage),
    stage= case_when(
      stage == "1" ~ "1a",
      stage == "2" ~ "1b",
      stage == "3" ~ "2a",
      stage == "4" ~ "2b",
      stage == "5" ~ "3a",
      stage == "6" ~ "3b",
      stage == "7" ~ "IV",
      TRUE ~ stage
    )
  ) %>% 
  mutate(
    histology = str_trim(histology),
    histology = case_when(
      histology == "1" ~ "SSC",
      histology == "2" ~ "unspecified_AC",
      histology == "3" ~ "Large cell",
      TRUE ~ histology
    )
  ) %>% 
  mutate(
    Smoking = str_trim(Smoking),
    Smoking = case_when(
      Smoking == "1" ~ "current",
      Smoking == "2" ~ "ex",
      Smoking == "3" ~ "never",
      TRUE ~ Smoking
    )
    
  )

#mutating the condition column into a readable format
metadata_new <- metadata_new %>%
  mutate(
    condition = case_when(
      str_ends(condition, "T") ~ "tumor",
      TRUE ~ "normal" 
    )
  )

#some additional pre-processing of metadata for correct alignment
metadata_new$samples <- gsub(".*_", "", metadata_new$samples)
#checking alignment
all(metadata_new$samples == colnames(count_raw))
all(colnames(count_raw) == rownames(metadata_new))
#sorting 
metadata_new$samples <- sort(metadata_new$samples)
colnames(count_raw) <- sort(colnames(count_raw))
#setting samples column as row names of metadata for deseq
metadata_new$GSM <- rownames(metadata_new)
rownames(metadata_new) <-  NULL
metadata_new$GSM <- rownames(metadata_new)
metadata_new <- metadata_new[, -which(names(metadata_new) == "GSM")]
rownames(metadata_new) <-  metadata_new$samples
metadata_new <- metadata_new[, -which(names(metadata_new)=="samples")]

##Merging
#merging for Visualization 
merged_data <- cbind(metadata_new, t(count_raw))
#converting categorical variables into factors
metadata_new$gender <- as.factor(metadata_new$gender)
metadata_new$Smoking <- as.factor(metadata_new$Smoking)
metadata_new$histology <- as.factor(metadata_new$histology)
metadata_new$stage <- as.factor(metadata_new$stage)
metadata_new$condition <- as.factor(metadata_new$condition)

#saving the count matrix
write.csv(count_raw, "count_matrix.csv", row.names = TRUE)
#saving the metadata
write.csv(metadata, "metadata_raw.csv", row.names = TRUE)
write.csv(metadata_new, "metadata_processed.csv", row.names = TRUE)

