library(readxl)
library(dplyr)
library(tibble)
library(phyloseq)

file_path <- "../data/dataset_MISO_project.xlsx"

otu_sheet <- "OTU_table"
tax_sheet <- "Taxonomy"
sample_sheet <- "Sample_data"

otu_mat <- read_excel(file_path, sheet = otu_sheet)
tax_mat <- read_excel(file_path, sheet = tax_sheet)
samples_df <- read_excel(file_path, sheet = sample_sheet)

otu_mat <- otu_mat %>% tibble::column_to_rownames("OTU_ID")
tax_mat <- tax_mat %>% tibble::column_to_rownames("Tax_ID")
samples_df <- samples_df %>% tibble::column_to_rownames("Sample_ID")

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

metagenomics <- phyloseq(OTU, TAX, samples)