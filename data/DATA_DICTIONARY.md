# ORCHiD Study - De-identified Sample Metadata

## File: orchid_metadata_deidentified.csv

**Description:** De-identified sample metadata for 132 ovarian cancer patients included in final microbiome analyses.

**Last updated:** January 2026

---

## Data Dictionary

### Sample Identifiers
- **SampleID**: Anonymized sample identifier (e.g., V100012)
  - Format: Letter prefix + numeric code
  - No linkage to personal identifiers

### Sequencing Quality Metrics
- **NumReads**: Total number of sequencing reads per sample after quality filtering
- **Shannon**: Shannon diversity index (alpha diversity metric)
- **Observed**: Number of observed ASVs (amplicon sequence variants)

### Demographic Variables (Original)
- **Race**: Self-reported race/ethnicity
  - Values: "White", "Black"
  - Note: Study limited to non-Hispanic White and non-Hispanic Black participants

- **age**: Age at diagnosis (years)
  - Continuous variable
  - Range: 30-79 years (per study eligibility)

### Demographic Variables (Analysis-Ready)
- **Race_clean**: Cleaned race variable for analysis
  - Values: "White", "Black"

- **age_clean**: Cleaned continuous age variable (same as 'age')
  - Used in continuous age models

- **age_cat_clean**: Categorized age for stratified analyses
  - Values: "<50", "≥50"
  - Used as proxy for menopausal status

### Clinical Variables (Original)
- **Stage**: Cancer stage at diagnosis
  - Values: "Early", "Late"

- **Histology_Group**: Histologic subtype
  - Original classifications

### Clinical Variables (Analysis-Ready)
- **Stage_clean**: Cleaned stage variable
  - Values: "Early" (FIGO I-II), "Late" (FIGO III-IV)

- **Histology_clean**: Cleaned histology variable for analysis
  - Values: 
    - "Type II epithelial" (high-grade serous and aggressive subtypes)
    - "Type I epithelial / Other" (type I epithelial, endometrioid, mucinous, clear cell, non-epithelial)

---

## Sample Summary

- **Total samples:** 132
- **Race distribution:**
  - Non-Hispanic White: 110 (83.3%)
  - Non-Hispanic Black: 22 (16.7%)
- **Age distribution:**
  - <50 years: 22 (16.7%)
  - ≥50 years: 110 (83.3%)
- **Stage distribution:**
  - Early (I-II): 41 (31.1%)
  - Late (III-IV): 91 (68.9%)

---

## Quality Filtering

Samples included in this dataset met the following criteria:
- ≥10,000 sequencing reads after quality filtering
- Passed DADA2 quality control
- No antibiotic use in 30 days prior to collection
- No vaginal product use in 30 days prior to collection
- Successful vaginal swab collection

Excluded samples:
- 4 negative controls (<1,000 reads)
- 3 samples with insufficient depth (<10,000 reads)

---

## Privacy Protection

**De-identification measures:**
- All personal identifiers removed
- Sample IDs are registry-generated codes with no linkage to names, dates, or locations
- Sensitive socioeconomic variables (income, employment, education, marital status) excluded
- Exact diagnosis dates not included
- Geographic information limited to state-level registry participation

**Permitted use:**
- Research purposes only
- Must cite original publication
- Cannot attempt re-identification

---

## Citation

If using this metadata, please cite:

Deveaux A, et al. (2026). Vaginal Microbiome Composition, Diversity and Dysbiosis: The ORCHiD Study. *Nature Communications Medicine* (in review).

---

## Contact

For questions about this dataset:
- **Corresponding Author:** April Deveaux, MD, PhD
- **Email:** april.deveaux@duke.edu
- **Affiliation:** Department of Population Health Sciences, Duke University School of Medicine

---

## Notes

- **Extraction batch information:** Not included in this file. Extraction batch showed no significant effect on community composition (R²=0.004, p=0.74) and was not included in final models.

- **Missing variables:** Some variables collected in the parent ORCHiD study (e.g., treatment details, exact diagnosis dates, socioeconomic variables) are not included to protect participant privacy.

- **Analysis variables:** The "_clean" suffix variables (Race_clean, age_clean, age_cat_clean, Stage_clean, Histology_clean) are the versions used in all statistical analyses reported in the manuscript.
