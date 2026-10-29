# ORCHiD Vaginal Microbiome Analyses

Analysis code for vaginal microbiome composition, dysbiosis, and perceived stress in ovarian cancer patients from the ORCHiD study.

This repository supports two analyses:

- **Vaginal microbiome composition and dysbiosis:** the original topic-modeling analysis of 132 participants, published in *Communications Medicine*.
- **Perceived stress and vaginal Lactobacillus:** the PSS-4 analysis of 123 participants, described in the manuscript *Perceived Stress and Vaginal Lactobacillus Species Among Patients With Ovarian Cancer in the ORCHiD Study*.

Each analysis has its own inputs, models, and software requirements. Use the instructions in the corresponding section below.

## Perceived stress and vaginal Lactobacillus analysis

### Analysis file

`ORCHiD_PSS4_Microbiome_Analysis_Publish.Rmd`

The R Markdown file contains the stress analysis and its manuscript comparisons. It starts from prepared analysis datasets. Sequencing processing, species assignment, microbiome measures, the PSS-4 total score, and the parent topic model are produced upstream.

### Study and models

- **Analytic sample:** 123 participants: 103 non-Hispanic White and 20 non-Hispanic Black.
- **Exposure:** PSS-4 total score, standardized using the full analytic sample; estimates are per 1-SD increase (3.43 points).
- **Species detection:** Firth logistic regression for *L. crispatus*, *L. iners*, *L. gasseri*, and *L. jensenii*, adjusted for race and age. Benjamini–Hochberg correction is applied across these four species tests.
- **Community and genus measures:** linear regression for Shannon diversity, Firth logistic regression for Lactobacillus dominance, and beta regression for Lactobacillus-dominated topic membership, adjusted for race, age, stage, and histologic subtype.
- **Race analyses:** age-adjusted interaction models and separate age-adjusted models within each racial group.
- **Sensitivity analyses:** expanded species adjustment, omission of histologic subtype, Spearman correlation among samples with detected *L. iners*, leave-one-out analysis, and checks for departure from linearity.

P values for community, genus, race, and sensitivity analyses are unadjusted.

### Required inputs

Place these three files in one directory and supply its location through the `data_dir` parameter:

| File | Contents used |
|---|---|
| `complete_stress_dysbiosis_dataset.csv` | Participant identifier, PSS-4 total score, age, race, stage, histologic subtype, Shannon diversity, and Lactobacillus dominance |
| `final_species_enhanced_dataset.csv` | Matching participant identifier, four species-detection indicators, and *L. iners* relative abundance |
| `topic_model_results.rds` | Parent topic-model output containing a `gamma` table with `document`, `topic`, and `gamma` columns |

The Rmd documents the required column names and coding. These prepared inputs are separate from the metadata supplied for the original dysbiosis analysis.

### Software

The Rmd requires `tidyverse`, `broom`, `logistf`, `betareg`, `patchwork`, and the R-recommended package `splines`. Rendering also requires `rmarkdown`, `knitr`, and Pandoc; RStudio supplies Pandoc.

The original stress-analysis HTML records R 4.6.1, `logistf` 1.26.1, and `betareg` 3.2-6. The Rmd records the environment used for each new run in `sessionInfo.txt`.

Install the CRAN packages before rendering:

```r
install.packages(c(
  "tidyverse", "broom", "logistf", "betareg", "patchwork",
  "rmarkdown", "knitr"
))
```

### Running the analysis

From the repository root, with the three prepared input files in `data/`:

```r
rmarkdown::render(
  "ORCHiD_PSS4_Microbiome_Analysis_Publish.Rmd",
  params = list(
    data_dir = "data",
    output_dir = "analysis_outputs",
    fail_on_mismatch = TRUE
  )
)
```

Replace `"data"` with the directory containing the prepared inputs when needed.

### Outputs

The rendered HTML contains the analysis and comparison results. The `analysis_outputs/` directory contains:

- `Table1_participant_characteristics.csv`
- `Table2_species_and_community_measures.csv`
- `Table3_race_interaction_and_stratified.csv`
- `SupplementaryTable_S1_sensitivity.csv`
- `analytic_event_counts.csv`
- `check_against_manuscript.csv`
- `Figure1.png`
- `Figure2.png`
- `sessionInfo.txt`

With `fail_on_mismatch = TRUE`, rendering stops if the final verification identifies missing, mismatched, or unexpected results. Inspect `check_against_manuscript.csv` for details.

### Data access

The prepared PSS-4 analysis datasets are not included in this code deposit. For inquiries about access to these inputs, contact April Deveaux at april.deveaux@duke.edu.

Raw sequencing reads are available through NCBI BioProject [PRJNA1414650](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1414650). Sequencing reads alone do not supply the survey and clinical variables required by the stress analysis.

## Original vaginal microbiome composition and dysbiosis analysis

The instructions below describe the original topic-modeling scripts and their inputs.

### Overview

This repository contains analysis code for characterizing vaginal microbiome composition in post-treatment ovarian cancer patients using 16S rRNA gene sequencing and Latent Dirichlet Allocation (LDA) topic modeling. The study identifies distinct microbial community patterns and their associations with demographic and clinical variables.

### Citation

Deveaux A, Osazuwa-Peters OL, Kim YJ, et al. (2026). Characterizing the composition and diversity of the vaginal microbiome in ovarian cancer. *Communications Medicine*. https://doi.org/10.1038/s43856-026-01771-8

### Repository Contents

```
orchid-microbiome-analysis/
├── README.md                              # This file
├── 01_topic_model_selection_and_stability.R    # K optimization & stability
├── 02_final_lda_model_k7.R                     # Final LDA model fitting
├── 03_beta_regression_topic_associations.R     # Differential abundance analysis
├── session_info.txt                            # R package versions
├── data/                                       # Input data directory
│   ├── README.md                               # Data availability information
│   └── orchid_metadata_deidentified.csv        # De-identified sample metadata
└── results/                                    # Output directory (created by scripts, not in repo)
```

**Note:** The processed phyloseq object (`ps_filtered.rds`) is not included in this repository due to size and privacy considerations, but is available upon request.

### Study Design

- **Population:** 132 ovarian cancer patients (post-treatment)
- **Sample Type:** Self-collected vaginal swabs
- **Sequencing:** 16S rRNA V1-V3 region, Illumina MiSeq
- **Analysis:** LDA topic modeling with beta regression for associations
- **Key Covariates:** Age (continuous), race, cancer stage, histologic subtype

### Data Availability

#### Raw Sequencing Data
Raw FASTQ files are available from the NCBI Sequence Read Archive:
- **BioProject:** PRJNA1414650

#### Processed Data
De-identified sample metadata and processed data tables are available in the manuscript supplementary materials and this repository.

#### Study Metadata
De-identified sample metadata is available in this repository:
- **File:** `data/orchid_metadata_deidentified.csv`
- **Contents:** 
  - Sample identifiers (anonymized)
  - Demographic variables (Race_clean, age_clean)
  - Clinical variables (Stage_clean, Histology_clean)
  - Sequencing batch information
  - Exclusion criteria flags (antibiotic use, vaginal product use)
- **Note:** All personal identifiers have been removed to protect participant privacy. This metadata file corresponds to the 132 samples included in final analyses after quality filtering.

#### Data Not Included in Repository
The following data files are **not included** in this public repository but are available as described:
- **Processed phyloseq object** (`ps_filtered.rds`): Available upon reasonable request due to file size and data privacy considerations
- **Raw sequencing data**: Available from NCBI SRA (BioProject PRJNA1414650)
- **Results files**: Generated locally when running analysis scripts

### Requirements

#### R Version
- R version 4.4.0 or higher

#### Required R Packages

##### Topic Modeling
- `topicmodels` (v0.2.17) - LDA implementation
- `tidytext` (v0.4.3) - Topic model interpretation
- `ldatuning` (v1.0.2) - Model selection metrics

##### Microbiome Analysis  
- `phyloseq` (v1.48.0) - Microbiome data handling
- `vegan` (v2.6.4) - Diversity metrics (if needed for other analyses)

##### Statistical Analysis
- `betareg` (v3.2-4) - Beta regression for proportions

##### Data Manipulation & Visualization
- `tidyverse` (v2.0.0) - Data wrangling
- `ggplot2` (v3.5.2) - Visualization
- `cowplot` (v1.2.0) - Multi-panel figures

See `session_info.txt` for complete package versions and dependencies.

### Installation

```r
# Install CRAN packages
install.packages(c("topicmodels", "tidytext", "ldatuning", 
                   "betareg", "tidyverse", "ggplot2", "cowplot"))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("phyloseq", "vegan"))
```

### Usage

#### 1. Data Preparation

**Option A: Use NCBI data**
Download raw sequencing data from NCBI SRA (BioProject PRJNA1414650) and process using DADA2 or your preferred pipeline. The scripts expect a quality-filtered phyloseq object with:
- Samples with >10,000 reads
- Taxa with >0.001% relative abundance across all samples
- Sample metadata including: `Race_clean`, `Stage_clean`, `Histology_clean`, `age_clean`

**Option B: Request processed data**
Contact the corresponding author for the processed phyloseq object (`ps_filtered.rds`).

**Required file location:**
Save your phyloseq object as `data/ps_filtered.rds` before running analysis scripts.

#### 2. Run Analysis Scripts in Order

##### Script 1: Topic Model Selection and Stability
```r
source("01_topic_model_selection_and_stability.R")
```

**Purpose:** Determines optimal number of topics (k) using multiple optimization metrics and tests model stability across random initializations.

**Outputs:**
- `results/model_fit_metrics.csv`
- `results/stability_k7_results.csv`
- `results/figures/model_fit_metrics.pdf`

**Time:** ~10-15 minutes

##### Script 2: Final LDA Model
```r
source("02_final_lda_model_k7.R")
```

**Purpose:** Fits final LDA model with k=7 topics using seed=243 for reproducibility.

**Outputs:**
- `results/lda_k7_model.rds`
- `results/beta_matrix_topic_genus_associations.csv`
- `results/gamma_matrix_sample_topic_proportions.csv`
- `results/topic_reference_guide.csv`
- `results/metadata_with_topic_proportions.csv`

**Time:** ~5 minutes

##### Script 3: Beta Regression Analysis
```r
source("03_beta_regression_topic_associations.R")
```

**Purpose:** Tests associations between topic proportions and clinical/demographic variables using beta regression.

**Outputs:**
- `results/beta_regression_comprehensive_results.csv`
- `results/beta_regression_[race/age/stage/histology]_results.csv`
- `results/beta_regression_manuscript_table.csv`

**Time:** <1 minute

#### 3. Interpret Results

Key results files for manuscript:
1. **Topic characterization:** `topic_reference_guide.csv`
2. **Topic associations:** `beta_regression_manuscript_table.csv`
3. **Model validation:** `stability_k7_results.csv`

### Key Analysis Parameters

#### Quality Filtering
- Minimum reads per sample: **10,000**
- Minimum taxa relative abundance: **0.001%** (1e-5)
- Final sample size: **132**

#### Topic Modeling (LDA)
- Method: **Gibbs sampling**
- Topic range tested: **k = 3-10**
- Final model: **k = 7**
- Random seed: **243**
- Gibbs iterations: **2000**
- Burn-in: **1000**
- Thinning: **100**

#### Stability Testing
- Number of random initializations: **20**
- Mean beta correlation: **0.883**

#### Beta Regression
- Link function: **logit**
- Formula: `Topic_proportion ~ Race_clean + Stage_clean + Histology_clean + age_clean`
- Multiple testing correction: **Benjamini-Hochberg FDR**
- Significance threshold: **FDR q < 0.05**

#### Variables
- **Age:** Continuous (years)
- **Race:** Non-Hispanic Black (reference) vs Non-Hispanic White
- **Stage:** Early (I-II, reference) vs Late (III-IV)
- **Histology:** Type I/non-epithelial (reference) vs Type II epithelial

### Reproducibility Notes

1. **Random Seeds:** All analyses use fixed random seeds for reproducibility
2. **Package Versions:** See `session_info.txt` for exact versions
3. **File Paths:** Scripts use relative paths; set working directory to repository root
4. **Computational Time:** Total runtime ~20 minutes on standard laptop

### Troubleshooting

#### Common Issues

**Error: "Cannot find phyloseq object"**
- Ensure your processed data is saved as `data/ps_filtered.rds`
- Check that sample metadata includes required variables

**Error: "Beta regression convergence issues"**
- This may occur if topic proportions are too extreme
- The epsilon transformation should handle this, but check your gamma values

**Error: "Package not found"**
- Install missing packages using installation commands above
- Some packages require Bioconductor

#### Getting Help

For questions about:
- **Analysis methods:** Contact april.deveaux@duke.edu
- **Code issues:** Open an issue in this repository
- **Data access:** See Data Availability section

### License

This code is released under the MIT License. See LICENSE file for details.

### Acknowledgments

This work was supported by:
- National Cancer Institute (K01CA255408 and R37CA233777)
- Duke University Population Health Sciences
- Duke Microbiome Core Facility
- State cancer registries (New York, Kentucky, California, North Carolina, Maryland, Georgia, Texas)

### References

1. **Topic Modeling Approach:** Blei DM, Ng AY, Jordan MI. (2003). Latent Dirichlet allocation. *J Mach Learn Res* 3:993-1022.

2. **Microbiome Topic Modeling:** Shafiei M, et al. (2015). BioMico: A supervised Bayesian model for inference of microbial community structure. *Microbiome* 3:8.

3. **Beta Regression:** Ferrari S, Cribari-Neto F. (2004). Beta regression for modelling rates and proportions. *J Appl Stat* 31:799-815.

4. **Phyloseq Package:** McMurdie PJ, Holmes S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLoS ONE* 8(4):e61217.

### Contact

**Corresponding Author:**  
April Deveaux, MD, PhD  
Department of Population Health Sciences  
Duke University School of Medicine  
Durham, NC 27701  
Email: april.deveaux@duke.edu

### Last Updated

October 29, 2026
