# Gallbladder incidence sensitivity analysis

Code-only R workflow for case-definition sensitivity, demographic standardization,
continuous-age trends, descriptive population decomposition and secondary
competing-risk analyses of gallbladder malignancy.

No manuscript, patient records, SEER export, aggregate result dataset, figure,
credential or author contact information is included. All inputs must be obtained
separately under the applicable SEER access and data-use terms.

## Requirements and execution

The analysis environment used R 4.6.1 and SEER*Stat 9.0.43.0 exports from the
November 2025 SEER Research Data submission (2000-2023). Required R packages:
splines, survival and cmprsk; optional plotting also requires ggplot2, patchwork,
scales and ragg. The runner checks packages but does not install software.

Run from the repository directory, supplying your own local paths:

```text
Rscript tests/test_methods.R
Rscript run_analysis.R OUTPUT_DIRECTORY RAW_DIRECTORY --no-figures
```

Omit `--no-figures` to generate TIFF statistical plots and local PNG previews.
Choose a new or empty output directory outside both the code and raw directories.
Inputs are read-only. Execution logs and generated results remain outside this
repository. Consult INPUT_SCHEMA.md for filenames and field conventions.

## Workflow

The runner executes the scripts in this order:

1. `08_redesign_core_analysis.R`: morphology-based counts, denominator checks,
   initial age-standardized trends, continuous-age models and coding composition.
2. `09_competing_risk_analysis.R`: outcome cohorts, cumulative incidence,
   Fine-Gray and cause-specific Cox models.
3. `11_revision_analysis.R`: unrounded directly standardized age-cutoff estimates,
   spline/window sensitivity, cohort flow and sequential outcome adjustment.
4. `12_time_varying_mortality.R`: age-by-log-time outcome associations.
5. `13_ethnicity_analysis.R`: common-population age-only/joint standardization,
   race-specific slopes and adjusted continuous-age models.
6. `14_reframe_analysis.R`: morphology audit, paired model contrasts,
   5,000-draw Poisson decomposition, sex trends and stage incidence reconstruction.
7. `15_robustness_analysis.R`: zero-cell counts/weights, 50,000-draw Poisson and
   Jeffreys Gamma posterior sensitivity, SEER classification-based sensitivity,
   and Fay-Feuer gamma stage intervals.
8. `10_manuscript_figures.R` (optional): statistical plotting from local results.

The numeric script names and some output keys are retained for compatibility.
Scripts 11 and 13 provide the directly standardized trend comparisons used in the
final analysis. Script 15 provides the current gamma stage intervals; the Wald
intervals produced by script 14 are intermediate outputs.

## Definitions and interpretation

- D1: all malignant C23.9 tumors. D2: numeric morphology range 8010-8576.
- D3: numeric range 8140-8389; it includes neuroendocrine/mixed entities and is not
  a complete or pure adenocarcinoma definition.
- D4: adenocarcinoma NOS, 8140/3. Legacy `strict` keys identify this code only.
- The additional malignant adenocarcinoma-with-variants cohort uses the exact
  SEER Rare Cancer Classification 15.1 / recode 60 code list. Its internal key
  `Glandular_subset` is a compatibility identifier for this variant category.

Age-group APC and annual change in the younger/older rate ratio are different
estimands. Legacy `differential_apc` fields mean rate-ratio annual change, not an
arithmetic difference in APC percentage points. The cutoff table's `dispersion`
field denotes WLS residual variance, not quasi-Poisson dispersion.

Decomposition separates absolute period-rate differences into composition and
within-stratum components. It is descriptive, not causal. Population denominators
are fixed. Jeffreys Gamma posterior intervals are credible intervals and differ
from plug-in Poisson bootstrap confidence intervals. The primary decomposition
uses 5,000 draws; robustness resimulations use 50,000. Seeds are fixed in code.
Zero-cell and posterior checks are not a claim of frequentist interval coverage.

SEER 8 overlaps SEER 17 and is a registry sensitivity analysis, not independent
validation. Rates are not reporting-delay adjusted. Sex models use rounded
exported rates and are not jointly ethnicity-standardized. Outcome analyses do
not establish treatment effects. No direct replication of another publication is
claimed.

The in-memory tests verify mathematical properties, not study-data reproduction.
Full reproduction additionally requires the exact authorized input exports.

## Method sources

- SEER Rare Cancer Classification, category 15.1:
  https://seer.cancer.gov/seerstat/variables/seer/raresiterecode/
- Symmetric rate decomposition: doi:10.1080/01621459.1955.10501299.
- Fay-Feuer gamma intervals: doi:10.1002/(SICI)1097-0258(19970415)16:7<791::AID-SIM500>3.0.CO;2-#.

No manuscript-building, author-management or submission-package scripts are
included. No automatic data download or network operation is performed by the R
workflow. No reuse license has been selected.
