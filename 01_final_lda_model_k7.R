################################################################################
# FINAL LDA MODEL: k=7 TOPICS
# 
# Purpose: Fit final LDA topic model with k=7 topics using seed=243
#
# This script implements the final topic modeling analysis reported in the 
# manuscript Methods section: "LDA was fitted to the data using the topicmodels
# package in R with Gibbs sampling method... The model was run with a fixed 
# seed (243) for reproducibility."
#
# Outputs:
#   - Beta matrix (topic-genus associations)
#   - Gamma matrix (sample-topic proportions)
#   - Topic characterizations and naming
#   - Phyloseq object with topic data for downstream analyses
#
# Author: April Deveaux
# Date: January 2026
# For: ORCHiD Vaginal Microbiome Study
################################################################################

# Required packages
library(topicmodels)
library(phyloseq)
library(tidyverse)
library(tidytext)

cat("================================================================================\n")
cat("FINAL LDA MODEL: k=7 TOPICS\n")
cat("================================================================================\n\n")

################################################################################
# STEP 1: LOAD AND PREPARE DATA
################################################################################

cat("Step 1: Loading and preparing phyloseq data...\n\n")

# Load quality-filtered phyloseq object
# This should be your phyloseq object after quality filtering:
#   - Samples with >10,000 reads
#   - Taxa with >0.001% relative abundance
ps <- readRDS("data/ps_filtered.rds")  # Update path to your data

# Quality filtering parameters (should already be applied, but documenting here)
# ps2 <- subset_samples(ps, NumReads > 10000)
# minTotRelAbun <- 1e-5  # 0.001% threshold
# x <- taxa_sums(ps2)
# keepTaxa <- (x / sum(x)) > minTotRelAbun
# ps2 <- prune_taxa(keepTaxa, ps2)

# Aggregate to genus level
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

# Create genus-named OTU table
taxa_names(ps_genus) <- tax_table(ps_genus)[, "Genus"]

cat("Data prepared:\n")
cat("  Samples:", nsamples(ps_genus), "\n")
cat("  Genera:", ntaxa(ps_genus), "\n\n")

################################################################################
# STEP 2: PREPARE COUNT MATRIX FOR LDA
################################################################################

cat("Step 2: Creating count matrix for LDA...\n\n")

# Extract count matrix (samples × genera)
# LDA requires integer counts, not relative abundances
count_matrix <- as.data.frame(otu_table(ps_genus))

# Verify format
cat("Count matrix dimensions:", nrow(count_matrix), "genera ×", ncol(count_matrix), "samples\n")
cat("Total counts:", sum(count_matrix), "\n\n")

################################################################################
# STEP 3: FIT FINAL LDA MODEL (k=7, seed=243)
################################################################################

cat("================================================================================\n")
cat("Step 3: Fitting LDA model with k=7 topics (seed=243)\n")
cat("================================================================================\n\n")

cat("Running LDA with Gibbs sampling...\n")
cat("This may take several minutes...\n\n")

# Fit LDA model with fixed seed for reproducibility
# Parameters match Methods section
set.seed(243)
lda_k7 <- LDA(count_matrix, 
              k = 7, 
              method = "Gibbs",
              control = list(
                seed = 243,
                iter = 2000,    # Number of Gibbs iterations
                burnin = 1000,  # Burn-in period
                thin = 100      # Thinning interval
              ))

cat("LDA model fitting complete!\n\n")

################################################################################
# STEP 4: EXTRACT BETA AND GAMMA MATRICES
################################################################################

cat("Step 4: Extracting beta and gamma matrices...\n\n")

# Beta matrix: Topic-term (topic-genus) distributions
# Shows which genera define each topic
# High beta = genus is highly associated with that topic
beta_df <- tidy(lda_k7, matrix = "beta")

cat("Beta matrix extracted:\n")
cat("  Dimensions:", nrow(beta_df), "topic-genus pairs\n")
cat("  Range:", round(min(beta_df$beta), 4), "-", round(max(beta_df$beta), 4), "\n\n")

# Gamma matrix: Document-topic (sample-topic) distributions  
# Shows topic mixture for each sample
# These are the proportions used in downstream beta regression
gamma_df <- tidy(lda_k7, matrix = "gamma") %>%
  arrange(document, topic)

cat("Gamma matrix extracted:\n")
cat("  Dimensions:", nrow(gamma_df), "sample-topic pairs\n")
cat("  Range:", round(min(gamma_df$gamma), 4), "-", round(max(gamma_df$gamma), 4), "\n\n")

################################################################################
# STEP 5: CHARACTERIZE TOPICS BY DOMINANT GENERA
################################################################################

cat("Step 5: Identifying dominant genera for each topic...\n\n")

# Find dominant genus (highest beta) for each topic
dominant_genera <- beta_df %>%
  group_by(topic) %>%
  slice_max(beta, n = 1) %>%
  ungroup() %>%
  arrange(topic)

cat("Dominant genera per topic:\n")
print(as.data.frame(dominant_genera))
cat("\n")

# Get top 3 genera per topic for more complete characterization
top3_genera <- beta_df %>%
  group_by(topic) %>%
  slice_max(beta, n = 3) %>%
  mutate(rank = row_number()) %>%
  arrange(topic, rank)

cat("Top 3 genera per topic:\n")
print(as.data.frame(top3_genera))
cat("\n")

################################################################################
# STEP 6: CREATE TOPIC REFERENCE GUIDE
################################################################################

cat("Step 6: Creating topic reference guide...\n\n")

# Create interpretive names based on dominant genera and known vaginal health states
# These names should match those reported in your manuscript
topic_reference <- data.frame(
  Topic_Number = 1:7,
  Short_Name = c(
    "Streptococcus-Mixed",
    "Proteus-Dominated", 
    "Diverse-Community",
    "Anaerobic-Dysbiosis",
    "Lactobacillus-Dominated",
    "Pseudomonas-Dominated",
    "BV-Type"
  ),
  Full_Name = c(
    "Streptococcus Mixed Community",
    "Proteus Dominated",
    "Diverse Intermediate Community",
    "Anaerobic Dysbiosis (Peptoniphilus)",
    "Lactobacillus Healthy",
    "Pseudomonas Environmental",
    "Bacterial Vaginosis Type (Prevotella)"
  ),
  Dominant_Genus = dominant_genera$term,
  Health_Status = c(
    "Moderate Dysbiosis",
    "Severe Dysbiosis",
    "Intermediate",
    "Transitional/Dysbiotic",
    "Healthy",
    "Pathogenic/Environmental",
    "BV-Associated Dysbiosis"
  )
)

cat("Topic Reference Guide:\n")
print(topic_reference)
cat("\n")

################################################################################
# STEP 7: CREATE PHYLOSEQ OBJECT WITH TOPIC DATA
################################################################################

cat("Step 7: Creating phyloseq object with topic proportions...\n\n")

# Reshape gamma to wide format (samples × topics)
gamma_wide <- gamma_df %>%
  pivot_wider(names_from = topic, 
              values_from = gamma, 
              names_prefix = "Topic_")

# Add sample metadata
sample_metadata <- data.frame(sample_data(ps_genus))
sample_metadata$SampleID <- rownames(sample_metadata)

# Merge topic proportions with metadata
metadata_with_topics <- sample_metadata %>%
  left_join(gamma_wide, by = c("SampleID" = "document"))

# Create new phyloseq object with topic proportions as "taxa"
# This allows using phyloseq functions for topic-level analyses

# Create OTU table from topic proportions
topic_counts <- gamma_wide %>%
  column_to_rownames("document") %>%
  t() %>%
  as.matrix()

# Create taxonomy table for topics
topic_tax <- data.frame(
  Kingdom = "Topic",
  Phylum = paste0("Topic", 1:7),
  Class = topic_reference$Short_Name,
  Order = topic_reference$Full_Name,
  Family = topic_reference$Dominant_Genus,
  Genus = topic_reference$Health_Status,
  row.names = rownames(topic_counts)
)

# Build phyloseq object
ps_topics <- phyloseq(
  otu_table(topic_counts, taxa_are_rows = TRUE),
  tax_table(as.matrix(topic_tax)),
  sample_data(metadata_with_topics %>% column_to_rownames("SampleID"))
)

cat("Phyloseq object with topic data created:\n")
print(ps_topics)
cat("\n")

################################################################################
# STEP 8: SAVE OUTPUTS
################################################################################

cat("Step 8: Saving outputs...\n\n")

# Create output directory
dir.create("results", showWarnings = FALSE)

# Save LDA model object
saveRDS(lda_k7, "results/lda_k7_model.rds")
cat("Saved: results/lda_k7_model.rds\n")

# Save beta matrix
write.csv(beta_df, "results/beta_matrix_topic_genus_associations.csv", row.names = FALSE)
cat("Saved: results/beta_matrix_topic_genus_associations.csv\n")

# Save gamma matrix
write.csv(gamma_df, "results/gamma_matrix_sample_topic_proportions.csv", row.names = FALSE)
cat("Saved: results/gamma_matrix_sample_topic_proportions.csv\n")

# Save topic reference guide
write.csv(topic_reference, "results/topic_reference_guide.csv", row.names = FALSE)
cat("Saved: results/topic_reference_guide.csv\n")

# Save dominant genera
write.csv(dominant_genera, "results/dominant_genera_per_topic.csv", row.names = FALSE)
cat("Saved: results/dominant_genera_per_topic.csv\n")

# Save top 3 genera
write.csv(top3_genera, "results/top3_genera_per_topic.csv", row.names = FALSE)
cat("Saved: results/top3_genera_per_topic.csv\n")

# Save phyloseq object with topics
saveRDS(ps_topics, "results/ps_topics.rds")
cat("Saved: results/ps_topics.rds\n")

# Save metadata with topic proportions
write.csv(metadata_with_topics, "results/metadata_with_topic_proportions.csv", row.names = FALSE)
cat("Saved: results/metadata_with_topic_proportions.csv\n\n")

################################################################################
# STEP 9: SUMMARY STATISTICS
################################################################################

cat("================================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("================================================================================\n\n")

# Sample dominance patterns
cat("Sample dominance by topic:\n")
dominant_topics <- gamma_df %>%
  group_by(document) %>%
  slice_max(gamma, n = 1) %>%
  ungroup()

topic_dominance_counts <- table(dominant_topics$topic)
cat("\n")
for(i in 1:7) {
  cat("  Topic", i, "(", topic_reference$Short_Name[i], "):",
      topic_dominance_counts[i], "samples (",
      round(topic_dominance_counts[i] / nsamples(ps_genus) * 100, 1), "%)\n")
}
cat("\n")

# Mean topic proportions across all samples
cat("Mean topic proportions across cohort:\n")
mean_proportions <- gamma_df %>%
  group_by(topic) %>%
  summarise(mean_gamma = mean(gamma),
            sd_gamma = sd(gamma),
            min_gamma = min(gamma),
            max_gamma = max(gamma))

print(as.data.frame(mean_proportions))
cat("\n")

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Final LDA model (k=7, seed=243) successfully fit and outputs saved.\n\n")

cat("Key outputs for downstream analyses:\n")
cat("  1. Gamma matrix (topic proportions) → Use for beta regression\n")
cat("  2. Beta matrix (topic-genus associations) → Use for topic interpretation\n")
cat("  3. Phyloseq object with topics → Use for visualization\n")
cat("  4. Topic reference guide → Use for manuscript tables\n\n")

cat("Next step: Run 03_beta_regression_topic_associations.R\n\n")

################################################################################
# END OF SCRIPT
################################################################################
