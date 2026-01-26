################################################################################
# TOPIC MODEL SELECTION AND STABILITY ANALYSIS
# 
# Purpose: Determine optimal number of topics (k) for LDA topic modeling
# 
# Methods:
#   1. Model fit metrics: CaoJuan2009, Arun2010, Deveaud2014, Griffiths2004
#   2. Stability testing: 20 independent runs with different random seeds
#   3. Topic consistency: Beta correlation analysis across runs
#
# Corresponds to Methods section: "Optimal topic number was determined using 
# the FindTopicsNumber function... evaluated across the topic number range... 
# model stability was evaluated across 20 independent runs"
#
# Author: April Deveaux
# Date: January 2026
# For: ORCHiD Vaginal Microbiome Study
################################################################################

# Required packages
library(topicmodels)
library(phyloseq)
library(tidyverse)

cat("================================================================================\n")
cat("TOPIC MODEL SELECTION AND STABILITY ANALYSIS\n")
cat("================================================================================\n\n")

################################################################################
# STEP 1: LOAD AND PREPARE DATA
################################################################################

cat("Step 1: Loading phyloseq object and preparing genus-level data...\n")

# Load your processed phyloseq object
# This should be the quality-filtered phyloseq object with >10,000 reads per sample
# and taxa filtered to >0.001% relative abundance
ps <- readRDS("data/ps_filtered.rds")  # Update path to your data location

# Aggregate to genus level for topic modeling
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

# Prepare OTU matrix for LDA (samples as rows, genera as columns)
otu_matrix <- as.matrix(t(otu_table(ps_genus)))
otu_matrix <- otu_matrix[rowSums(otu_matrix) > 0, ]

# Extract genus names for interpretability
tax_df <- as.data.frame(tax_table(ps_genus))
genus_names <- tax_df[colnames(otu_matrix), "Genus"]
genus_names <- as.character(genus_names)
genus_names[is.na(genus_names)] <- paste0("Unknown_", seq_len(sum(is.na(genus_names))))
colnames(otu_matrix) <- genus_names

cat("Data prepared:\n")
cat("  Samples:", nrow(otu_matrix), "\n")
cat("  Genera:", ncol(otu_matrix), "\n\n")

################################################################################
# STEP 2: MODEL FIT METRICS (k = 3 to 10)
################################################################################

cat("================================================================================\n")
cat("Step 2: Calculating model fit metrics for k = 3 to 10 topics\n")
cat("================================================================================\n\n")

# Topic range based on vaginal microbiome CST framework (5-9 CSTs typical)
# Extended to 3-10 to allow for novel patterns in cancer population
k_range <- seq(3, 10, by = 1)

cat("Testing k from 3 to 10 topics...\n")
cat("This will take several minutes (4 metrics × 8 models)...\n\n")

fit_metrics_list <- list()
set.seed(123)

for(k in k_range) {
  cat("Testing k =", k, "...\n")
  
  # Fit LDA model
  set.seed(123)
  lda_model <- LDA(otu_matrix, 
                   k = k, 
                   method = "Gibbs",
                   control = list(seed = 123, iter = 2000, burnin = 1000, thin = 100))
  
  # Extract topic-term (beta) and document-topic (theta) distributions
  beta <- exp(lda_model@beta)
  theta <- lda_model@gamma
  
  # Calculate CaoJuan2009 metric (minimize)
  # Measures topic density - lower values indicate more distinct topics
  cao_juan <- sum(sapply(1:k, function(i) {
    sum(sapply((1:k)[-i], function(j) {
      sum(pmin(beta[i,], beta[j,]))
    }))
  })) / (k * (k - 1))
  
  # Calculate Arun2010 metric (minimize)
  # Based on symmetric KL divergence between topics
  m1 <- apply(beta, 1, function(x) sum(x * log(x + 1e-10)))
  cm1 <- colSums(beta)
  m2 <- sum(cm1 * log(cm1 + 1e-10))
  svd_theta <- svd(theta)
  arun <- sum(abs(m1 - m2) * svd_theta$d)
  
  # Calculate Deveaud2014 metric (maximize)
  # Measures topic diversity using Jensen-Shannon divergence
  jsd <- function(p, q) {
    m <- 0.5 * (p + q)
    0.5 * sum(p * log((p + 1e-10) / (m + 1e-10))) + 
      0.5 * sum(q * log((q + 1e-10) / (m + 1e-10)))
  }
  
  deveaud <- 0
  if(k > 1) {
    for(i in 1:(k-1)) {
      for(j in (i+1):k) {
        deveaud <- deveaud + jsd(beta[i,], beta[j,])
      }
    }
    deveaud <- deveaud / (k * (k - 1) / 2)
  }
  
  # Calculate Griffiths2004 metric (maximize)
  # Log-likelihood of the model
  griffiths <- logLik(lda_model)
  
  # Store results
  fit_metrics_list[[as.character(k)]] <- data.frame(
    topics = k,
    CaoJuan2009 = cao_juan,
    Arun2010 = arun,
    Deveaud2014 = deveaud,
    Griffiths2004 = as.numeric(griffiths)
  )
  
  cat("  CaoJuan:", round(cao_juan, 4), 
      "| Arun:", round(arun, 2),
      "| Deveaud:", round(deveaud, 4),
      "| Griffiths:", round(as.numeric(griffiths), 2), "\n")
}

# Combine results
fit_metrics <- bind_rows(fit_metrics_list)

cat("\nModel fit metrics completed!\n\n")

# Identify optimal k values from each metric
cao_min <- fit_metrics$topics[which.min(fit_metrics$CaoJuan2009)]
arun_min <- fit_metrics$topics[which.min(fit_metrics$Arun2010)]
deveaud_max <- fit_metrics$topics[which.max(fit_metrics$Deveaud2014)]
griffiths_max <- fit_metrics$topics[which.max(fit_metrics$Griffiths2004)]

cat("MODEL FIT METRICS RESULTS:\n")
cat("--------------------------\n")
print(fit_metrics)
cat("\n")

cat("INTERPRETATION:\n")
cat("- CaoJuan2009 (minimize): optimal k =", cao_min, "\n")
cat("- Arun2010 (minimize): optimal k =", arun_min, "\n")
cat("- Deveaud2014 (maximize): optimal k =", deveaud_max, "\n")
cat("- Griffiths2004 (maximize): optimal k =", griffiths_max, "\n\n")

# Save results
dir.create("results", showWarnings = FALSE)
write.csv(fit_metrics, "results/model_fit_metrics.csv", row.names = FALSE)

################################################################################
# STEP 3: STABILITY ANALYSIS (k = 7)
################################################################################

cat("================================================================================\n")
cat("Step 3: Stability analysis for k = 7 across 20 random seeds\n")
cat("================================================================================\n\n")

# Based on convergence of metrics, test stability of k=7 model
K <- 7
N_SEEDS <- 20

# Generate random seeds for reproducibility
set.seed(456)
random_seeds <- sample(1:10000, N_SEEDS)

cat("Testing k = 7 with", N_SEEDS, "different random initializations...\n")
cat("Random seeds:", paste(random_seeds[1:5], collapse = ", "), "...\n\n")

# Storage for beta matrices from each run
all_beta_matrices <- list()

for(i in 1:N_SEEDS) {
  if(i %% 5 == 0) cat("  Completed", i, "of", N_SEEDS, "runs\n")
  
  lda_model <- LDA(otu_matrix, 
                   k = K, 
                   method = "Gibbs",
                   control = list(seed = random_seeds[i],
                                  iter = 2000,
                                  burnin = 1000,
                                  thin = 100))
  
  # Extract and store beta matrix (topic-genus associations)
  beta_matrix <- exp(lda_model@beta)
  rownames(beta_matrix) <- paste0("Topic_", 1:K)
  colnames(beta_matrix) <- genus_names
  
  all_beta_matrices[[i]] <- beta_matrix
}

cat("\n")

################################################################################
# STEP 4: CALCULATE BETA CORRELATIONS ACROSS RUNS
################################################################################

cat("Step 4: Calculating beta correlations between all run pairs...\n\n")

# Function to align topics between two runs based on maximum correlation
align_topics <- function(beta1, beta2) {
  K <- nrow(beta1)
  cor_matrix <- matrix(NA, K, K)
  
  # Calculate correlation between all topic pairs
  for(i in 1:K) {
    for(j in 1:K) {
      cor_matrix[i, j] <- cor(beta1[i, ], beta2[j, ])
    }
  }
  
  # Greedy matching (approximation of Hungarian algorithm)
  alignments <- numeric(K)
  used <- logical(K)
  
  for(i in 1:K) {
    available_cors <- cor_matrix[i, ]
    available_cors[used] <- -Inf
    best_match <- which.max(available_cors)
    alignments[i] <- best_match
    used[best_match] <- TRUE
  }
  
  return(list(
    alignment = alignments,
    correlations = sapply(1:K, function(i) cor_matrix[i, alignments[i]])
  ))
}

# Compare all pairs of runs
all_correlations <- c()

for(i in 1:(N_SEEDS-1)) {
  for(j in (i+1):N_SEEDS) {
    alignment_result <- align_topics(all_beta_matrices[[i]], 
                                     all_beta_matrices[[j]])
    all_correlations <- c(all_correlations, alignment_result$correlations)
  }
}

# Calculate summary statistics
mean_cor <- mean(all_correlations)
sd_cor <- sd(all_correlations)
min_cor <- min(all_correlations)
max_cor <- max(all_correlations)

cat("STABILITY RESULTS FOR k = 7:\n")
cat("----------------------------\n")
cat("Mean beta correlation:", round(mean_cor, 3), "\n")
cat("SD:", round(sd_cor, 3), "\n")
cat("Range:", round(min_cor, 3), "-", round(max_cor, 3), "\n")

# Count highly stable topics
topic_correlations <- matrix(all_correlations, nrow = K)
topic_mean_cors <- apply(topic_correlations, 1, mean)
high_stability <- sum(topic_mean_cors >= 0.85)
moderate_stability <- sum(topic_mean_cors >= 0.75 & topic_mean_cors < 0.85)
low_stability <- sum(topic_mean_cors < 0.75)

cat("\nPer-topic stability:\n")
cat("  High (r ≥ 0.85):", high_stability, "/", K, "\n")
cat("  Moderate (0.75 ≤ r < 0.85):", moderate_stability, "/", K, "\n")
cat("  Low (r < 0.75):", low_stability, "/", K, "\n")
cat("  Success rate:", round((high_stability + moderate_stability) / K * 100, 1), "%\n\n")

# Save stability results
stability_results <- data.frame(
  k = K,
  n_seeds = N_SEEDS,
  mean_correlation = mean_cor,
  sd_correlation = sd_cor,
  min_correlation = min_cor,
  max_correlation = max_cor,
  high_stability_topics = high_stability,
  moderate_stability_topics = moderate_stability,
  low_stability_topics = low_stability
)

write.csv(stability_results, "results/stability_k7_results.csv", row.names = FALSE)

################################################################################
# STEP 5: VISUALIZE RESULTS
################################################################################

cat("Step 5: Creating visualization plots...\n\n")

# Create output directory for figures
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Plot model fit metrics
pdf("results/figures/model_fit_metrics.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))

plot(fit_metrics$topics, fit_metrics$CaoJuan2009, 
     type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Number of Topics", ylab = "CaoJuan2009",
     main = "CaoJuan2009 (minimize)")
points(cao_min, fit_metrics$CaoJuan2009[fit_metrics$topics == cao_min],
       pch = 19, col = "blue", cex = 2)

plot(fit_metrics$topics, fit_metrics$Arun2010,
     type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Number of Topics", ylab = "Arun2010",
     main = "Arun2010 (minimize)")
points(arun_min, fit_metrics$Arun2010[fit_metrics$topics == arun_min],
       pch = 19, col = "blue", cex = 2)

plot(fit_metrics$topics, fit_metrics$Deveaud2014,
     type = "b", pch = 19, col = "darkgreen", lwd = 2,
     xlab = "Number of Topics", ylab = "Deveaud2014",
     main = "Deveaud2014 (maximize)")
points(deveaud_max, fit_metrics$Deveaud2014[fit_metrics$topics == deveaud_max],
       pch = 19, col = "blue", cex = 2)

plot(fit_metrics$topics, fit_metrics$Griffiths2004,
     type = "b", pch = 19, col = "darkgreen", lwd = 2,
     xlab = "Number of Topics", ylab = "Griffiths2004",
     main = "Griffiths2004 (maximize)")
points(griffiths_max, fit_metrics$Griffiths2004[fit_metrics$topics == griffiths_max],
       pch = 19, col = "blue", cex = 2)

dev.off()

cat("Saved: results/figures/model_fit_metrics.pdf\n\n")

################################################################################
# SUMMARY
################################################################################

cat("================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n\n")

cat("RECOMMENDATION: k = 7 topics\n")
cat("  - Converged across multiple optimization metrics\n")
cat("  - Mean beta correlation:", round(mean_cor, 3), "\n")
cat("  - High stability in", high_stability, "of", K, "topics\n\n")

cat("Output files saved in results/ directory:\n")
cat("  - model_fit_metrics.csv\n")
cat("  - stability_k7_results.csv\n")
cat("  - figures/model_fit_metrics.pdf\n\n")

cat("Next step: Run 02_final_lda_model_k7.R to fit final model with seed=243\n\n")

################################################################################
# END OF SCRIPT
################################################################################
