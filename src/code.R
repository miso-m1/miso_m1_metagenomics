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

# metagenomics <- rarefy_even_depth(metagenomics, sample.size = min(sample_sums(metagenomics)))
# metagenomics

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

# function to test if rarefaction is needed
# rarefaction is only needed if the sample sizes are different, so we can check the sample sizes first
# samples size ~ 100
# The range of samples sizes should be checked to see if they are different enough to warrant rarefaction. If the sample sizes are similar, rarefaction may not be necessary.
# a good range is if the sample sizes differ by more than 10-20%. If the sample sizes are within this range, rarefaction may not be necessary. However, if the sample sizes differ significantly (e.g., one sample has 100 reads and another has 1000 reads), rarefaction may be needed to ensure that the samples are comparable.
rarefaction_needed <- function(physeq) {
  sample_sizes <- sample_sums(physeq)
  size_range <- range(sample_sizes)
  size_diff <- size_range[2] - size_range[1]
  size_diff_percent <- (size_diff / size_range[1]) * 100
  return(size_diff_percent > 5 && size_diff_percent < 5) # if the difference is greater than and less than 20%, rarefaction is needed
}

# return the sample sizes that are needed for rarefaction
if (rarefaction_needed(metagenomics)) {
  print("Rarefaction is needed")
  min_sample_size <- min(sample_sums(metagenomics))
  metagenomics <- rarefy_even_depth(metagenomics, sample.size = min_sample_size)
} else {
  print("Rarefaction is not needed")
}




