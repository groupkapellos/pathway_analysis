# Load libraries
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(tidyr)

setwd("C:/Users/LisaNicoleBörner/Downloads/")

# Step 1: Read Excel File
file_path <- "pathway_gene_file.xlsx"
df <- read_excel(file_path)

# Step 2: Function to clean and split gene strings into unique vectors
clean_split <- function(cell) {
  if (is.na(cell)) return(character(0))
  genes <- str_split(cell, ",\\s*")[[1]]
  unique(trimws(genes))
}

# Step 3: Clean the data (remove duplicate genes within each cell)
df_cleaned <- df %>%
  mutate(
    Positive_Genes = sapply(Positive_Genes, function(x) paste(unique(clean_split(x)), collapse = ", ")),
    Negative_Genes = sapply(Negative_Genes, function(x) paste(unique(clean_split(x)), collapse = ", "))
  )

# Step 4: Create overlap table: genes found in both Positive and Negative lists
overlap_df <- df_cleaned %>%
  rowwise() %>%
  mutate(
    Positive_List = list(clean_split(Positive_Genes)),
    Negative_List = list(clean_split(Negative_Genes)),
    Overlapping_Genes = list(intersect(Positive_List, Negative_List)),
    Overlapping_Genes_Str = paste(sort(Overlapping_Genes), collapse = ", ")
  ) %>%
  ungroup() %>%
  select(pathway, Overlapping_Genes_Str)

# Step 5: Write both outputs to Excel
write_xlsx(df_cleaned, "cleaned_pathway_gene_file.xlsx")
write_xlsx(overlap_df, "overlapping_genes_per_pathway.xlsx")

