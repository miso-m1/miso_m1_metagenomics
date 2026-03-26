# ---
# title: "Statistical Analysis of Metagenomic Data"
# format: html
# editor: visual
# ---

library("phyloseq")
library("ggplot2")      # graphics
library("readxl")       # necessary to import the data from Excel file
library("dplyr")        # filter and reformat data frames
library("tibble")       # Needed for converting column to row names
library("readxl") # necessary to import the data from Excel file

otu_mat<- read_excel("data/dataset_MISO_project.xlsx", sheet = "OTU_table")
tax_mat<- read_excel("data/dataset_MISO_project.xlsx", sheet = "Taxonomy")
samples_df <- read_excel("data/dataset_MISO_project.xlsx", sheet = "Sample_data")

# Convert data frames to matrices and set row names. Phyloseq objects need to have row.names
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU_ID")
tax_mat <- tax_mat %>%
  tibble::column_to_rownames("Tax_ID")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("Sample_ID")

#Transform into matrixes otu and tax tables (sample table can be left as data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

# Transform to phyloseq objects
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

metagenomics <- phyloseq(OTU, TAX, samples)
metagenomics

sample_names(metagenomics)
rank_names(metagenomics)
sample_variables(metagenomics)

metagenomics <- rarefy_even_depth(metagenomics, sample.size = min(sample_sums(metagenomics)))
metagenomics

# diseases bar plot
diseases <- samples_df$disease
country <- samples_df$country

# diseases by country count
disease_country_count <- table(diseases, country)
disease_country_count

# sample count
# country_count <- table(country)
# disease_count <- table(diseases)

# bar plot of country count
barplot(country_count, main = "Sample count by country", xlab = "Country", ylab = "Count", col = "blue")
country_count

# # scatterplot of disease by country
# plot(diseases, country, main = "Disease by Country", xlab = "country", ylab = "disease", col = "red")

# ANOVA test to check if basteria abundance differs by disease status
# Extract the abundance of a specific bacteria (e.g., "Bacteroides fragilis abundance")







