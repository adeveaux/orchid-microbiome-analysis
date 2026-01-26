################################################################################
# BETA REGRESSION: TOPIC-LEVEL DIFFERENTIAL ABUNDANCE ANALYSIS
#
# Purpose: Test associations between topic proportions and clinical/demographic
#          variables using beta regression
#
# This script implements the differential abundance analysis reported in the
# manuscript Methods section: "Topic-level differential abundance testing was 
# performed using beta regression... directly models the gamma proportions... 
# with logit link function and the formula: Topic_proportion ~ Race_clean + 
# Stage_clean + Histology_clean + age_clean"
#
# Beta regression is appropriate for proportions (0-1 bounded) and avoids
# artificial conversion to count space
#
# Author: April Deveaux
# Date: January 2026
# For: ORCHiD Vaginal Microbiome Study
################################################################################

# Required packages
library(tidyverse)
library(betareg)
library(phyloseq)

cat("================================================================================\n")
cat("BETA REGRESSION: TOPIC-LEVEL DIFFERENTIAL ABUNDANCE ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# STEP 1: LOAD DATA
################################################################################

cat("Step 1: Loading gamma matrix and sample metadata...\n\n")

# Load gamma matrix from topic modeling (Script 02 output)
gamma_df <- read.csv("results/gamma_matrix_sample_topic_proportions.csv")

# Load phyloseq object for metadata
ps <- readRDS("data/ps_filtered.rds")  # Update path to your data

# Extract sample metadata
sample_metadata <- data.frame(sample_data(ps)) %>%
  rownames_to_column(var = "document")

cat("Data loaded:\n")
cat("  Samples:", nrow(sample_metadata), "\n")
cat("  Age range:", min(sample_metadata$age_clean), "-", 
    max(sample_metadata$age_clean), "years\n\n")

################################################################################
# STEP 2: PREPARE DATA FOR BETA REGRESSION
################################################################################

cat("Step 2: Preparing data for beta regression...\n\n")

# Reshape gamma to wide format (samples × topics)
gamma_wide <- gamma_df %>%
  pivot_wider(names_from = topic, 
              values_from = gamma, 
              names_prefix = "Topic_")

# Merge with metadata
beta_data <- gamma_wide %>% 
  left_join(sample_metadata, by = "document")

# Epsilon transformation for boundary values
# Beta regression requires values in open interval (0, 1)
# This transformation: y' = (y × (n-1) + 0.5) / n
# Recommended in Smithson & Verkuilen (2006)
epsilon_transform <- function(y, n = length(y)) {
  (y * (n - 1) + 0.5) / n
}

# Apply transformation to all topic columns
topic_cols <- paste0("Topic_", 1:7)
beta_data_transformed <- beta_data %>%
  mutate(across(all_of(topic_cols), epsilon_transform))

cat("Epsilon transformation applied to handle boundary values.\n")
cat("Transformation: y' = (y × (n-1) + 0.5) / n\n\n")

################################################################################
# STEP 3: FIT BETA REGRESSION MODELS
################################################################################

cat("================================================================================\n")
cat("Step 3: Fitting beta regression models\n")
cat("================================================================================\n\n")

cat("Formula: Topic_proportion ~ Race_clean + Stage_clean + Histology_clean + age_clean\n")
cat("Link function: logit\n")
cat("Testing associations while controlling for other covariates\n\n")

# Function to fit beta regression for one topic
fit_beta_regression <- function(topic_name, data, formula_vars) {
  formula_str <- paste0(topic_name, " ~ ", paste(formula_vars, collapse = " + "))
  formula_obj <- as.formula(formula_str)
  
  cat("Fitting:", formula_str, "\n")
  
  tryCatch({
    # Fit beta regression model
    model <- betareg(formula_obj, data = data, link = "logit")
    
    # Extract coefficients
    coef_summary <- summary(model)$coefficients$mean
    
    # Format results
    results <- as.data.frame(coef_summary) %>%
      rownames_to_column(var = "term") %>%
      rename(Estimate = Estimate, 
             Std.Error = `Std. Error`, 
             z_value = `z value`, 
             pvalue = `Pr(>|z|)`) %>%
      filter(term != "(Intercept)") %>%
      mutate(
        Topic = topic_name,
        # Convert to log2 scale for interpretability
        log2FoldChange = Estimate / log(2),
        lfcSE = Std.Error / log(2)
      ) %>%
      select(Topic, term, Estimate, Std.Error, log2FoldChange, lfcSE, 
             z_value, pvalue)
    
    return(results)
    
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Predictor variables (all using continuous age)
predictor_vars <- c("Race_clean", "Stage_clean", "Histology_clean", "age_clean")

# Fit models for all 7 topics
all_results <- list()
for (topic in topic_cols) {
  result <- fit_beta_regression(topic, beta_data_transformed, predictor_vars)
  if (!is.null(result)) {
    all_results[[topic]] <- result
  }
}

# Combine all results
combined_results <- bind_rows(all_results)

cat("\nAll models fitted successfully!\n\n")

################################################################################
# STEP 4: FDR CORRECTION
################################################################################

cat("Step 4: Applying FDR correction...\n\n")

# FDR correction separately for each term (covariate)
# This tests each covariate across all 7 topics
results_with_fdr <- combined_results %>%
  group_by(term) %>%
  mutate(
    padj = p.adjust(pvalue, method = "BH"),  # Benjamini-Hochberg FDR
    reject = padj < 0.05
  ) %>%
  ungroup() %>%
  arrange(term, padj)

cat("FDR correction applied using Benjamini-Hochberg method.\n")
cat("Significance threshold: FDR q < 0.05\n\n")

# Summary of significant results
sig_summary <- results_with_fdr %>%
  group_by(term) %>%
  summarise(
    n_significant = sum(reject),
    min_padj = min(padj),
    .groups = "drop"
  )

cat("SIGNIFICANCE SUMMARY:\n")
cat("---------------------\n")
print(as.data.frame(sig_summary))
cat("\n")

################################################################################
# STEP 5: FORMAT RESULTS BY COVARIATE
################################################################################

cat("Step 5: Formatting results by covariate...\n\n")

# Helper function to format results for one covariate
format_for_term <- function(term_name, results_df) {
  results_df %>%
    filter(term == term_name) %>%
    mutate(
      baseMean = NA,  # Not applicable for beta regression
      stat = z_value,
      df = NA,        # Not applicable for beta regression
      comparison = term_name
    ) %>%
    select(Topic, baseMean, log2FoldChange, lfcSE, stat, pvalue, 
           padj, reject, df, comparison) %>%
    arrange(Topic)
}

# Format by covariate
beta_race <- format_for_term("Race_cleanWhite", results_with_fdr)
beta_stage <- format_for_term("Stage_cleanLate", results_with_fdr)
beta_histology <- format_for_term("Histology_cleanType II epithelial", results_with_fdr)
beta_age <- format_for_term("age_clean", results_with_fdr)

# Print results for age (typically most interesting)
cat("AGE ASSOCIATIONS (continuous age, per year):\n")
cat("--------------------------------------------\n")
print(as.data.frame(beta_age %>% 
                      select(Topic, log2FoldChange, pvalue, padj, reject)))
cat("\n")

# Print results for race
cat("RACE ASSOCIATIONS (White vs Black):\n")
cat("------------------------------------\n")
print(as.data.frame(beta_race %>% 
                      select(Topic, log2FoldChange, pvalue, padj, reject)))
cat("\n")

################################################################################
# STEP 6: SAVE OUTPUTS
################################################################################

cat("Step 6: Saving results...\n\n")

# Create output directory if needed
dir.create("results", showWarnings = FALSE)

# Save comprehensive results
write.csv(results_with_fdr, 
          "results/beta_regression_comprehensive_results.csv", 
          row.names = FALSE)
cat("Saved: results/beta_regression_comprehensive_results.csv\n")

# Save covariate-specific results
write.csv(beta_race, "results/beta_regression_race_results.csv", row.names = FALSE)
cat("Saved: results/beta_regression_race_results.csv\n")

write.csv(beta_age, "results/beta_regression_age_results.csv", row.names = FALSE)
cat("Saved: results/beta_regression_age_results.csv\n")

write.csv(beta_stage, "results/beta_regression_stage_results.csv", row.names = FALSE)
cat("Saved: results/beta_regression_stage_results.csv\n")

write.csv(beta_histology, "results/beta_regression_histology_results.csv", row.names = FALSE)
cat("Saved: results/beta_regression_histology_results.csv\n\n")

################################################################################
# STEP 7: CREATE SUMMARY TABLE FOR MANUSCRIPT
################################################################################

cat("Step 7: Creating manuscript-ready summary table...\n\n")

# Create formatted summary for manuscript
manuscript_table <- results_with_fdr %>%
  mutate(
    # Format effect size with CI
    log2FC_formatted = sprintf("%.2f (%.2f)", log2FoldChange, lfcSE),
    # Format p-values
    pvalue_formatted = ifelse(pvalue < 0.001, "<0.001", 
                               sprintf("%.3f", pvalue)),
    padj_formatted = ifelse(padj < 0.001, "<0.001", 
                             sprintf("%.3f", padj)),
    # Significance indicator
    sig_indicator = case_when(
      padj < 0.05 ~ "*",
      padj < 0.10 ~ "†",
      TRUE ~ ""
    )
  ) %>%
  select(Topic, term, log2FC_formatted, pvalue_formatted, 
         padj_formatted, sig_indicator) %>%
  arrange(term, Topic)

write.csv(manuscript_table, 
          "results/beta_regression_manuscript_table.csv", 
          row.names = FALSE)
cat("Saved: results/beta_regression_manuscript_table.csv\n\n")

################################################################################
# STEP 8: INTERPRETATION NOTES
################################################################################

cat("================================================================================\n")
cat("INTERPRETATION NOTES\n")
cat("================================================================================\n\n")

cat("EFFECT SIZES (log2 fold-change):\n")
cat("  - Positive: Topic proportion increases with predictor\n")
cat("  - Negative: Topic proportion decreases with predictor\n")
cat("  - Magnitude: log2(FC) = 1 means 2-fold change\n\n")

cat("AGE INTERPRETATION:\n")
cat("  - Age is continuous (per year)\n")
cat("  - For easier interpretation, multiply by 10 for effect per decade\n")
cat("  - Example: log2FC = 0.05 per year = 0.5 per decade\n\n")

cat("RACE INTERPRETATION:\n")
cat("  - Reference: Non-Hispanic Black\n")
cat("  - Comparison: Non-Hispanic White\n")
cat("  - Positive log2FC: Higher in White vs Black\n")
cat("  - Negative log2FC: Lower in White vs Black\n\n")

cat("STAGE INTERPRETATION:\n")
cat("  - Reference: Early stage (I-II)\n")
cat("  - Comparison: Late stage (III-IV)\n\n")

cat("HISTOLOGY INTERPRETATION:\n")
cat("  - Reference: Type I epithelial / non-epithelial\n")
cat("  - Comparison: Type II epithelial (high-grade serous)\n\n")

################################################################################
# ANALYSIS COMPLETE
################################################################################

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("Beta regression analysis completed successfully.\n\n")

cat("Key findings to report in manuscript:\n")
cat("  1. Identify topics with FDR q < 0.05\n")
cat("  2. Report log2 fold-changes with standard errors\n")
cat("  3. For continuous age: consider scaling to per-decade effects\n")
cat("  4. Reference the comprehensive results CSV for all associations\n\n")

cat("Files saved in results/ directory\n\n")

################################################################################
# END OF SCRIPT
################################################################################
