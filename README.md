# Gallbladder cancer population-composition analysis

R code for age-only and joint age–ethnicity standardization of gallbladder cancer incidence in US SEER registries, paired APC comparisons, descriptive rate decomposition, and sex-stratified sensitivity analysis.

Software release: 1.2.0. Archived code: [10.5281/zenodo.23812431](https://doi.org/10.5281/zenodo.23812431).

This is a code-only repository. It contains no manuscript, source export, case record, population dataset, analytical result dataset, figure, or author contact details. Obtain the authorized inputs separately under the applicable SEER access and data-use terms. Do not place inputs or generated outputs inside this repository.

## Run

The analysis environment used R 4.6.1 and SEER*Stat 9.0.43.0. Required R packages are base/recommended `splines` and `MASS`. Optional figure generation uses `ggplot2` and `patchwork`. The workflow checks dependencies but installs nothing and performs no network requests.

```text
Rscript tests/test_methods.R
Rscript run_analysis.R OUTPUT_DIRECTORY RAW_DIRECTORY
Rscript run_analysis.R ANOTHER_OUTPUT_DIRECTORY RAW_DIRECTORY --figures
```

Use a new or empty output directory, separate from both the code and raw-data directories. The parent output directory must already exist. Each run reads three exports described in [INPUT_SCHEMA.md](INPUT_SCHEMA.md), reconstructs the analysis inputs, and writes private inputs, results and execution logs outside the repository. No precomputed derived dataset is required.

## Workflow

| Script | Purpose |
|---|---|
| `scripts/prepare_inputs.R` | Validate C2/P2 counts, retain zero-count cells, apply morphology definitions, calculate reference shares and WLS intervals, and assemble period cells. |
| `scripts/paired_trends.R` | Fit window-specific Poisson generators; calculate paired APC contrasts, reference-weight sensitivity, calendar diagnostics, gamma rate intervals and descriptive decomposition. |
| `scripts/uncertainty_checks.R` | Fit window-specific negative binomial generators and calculate annual-rate Fay–Feuer intervals. |
| `scripts/sex_analysis.R` | Validate P3 against P2 and compute sex-specific APCs, paired contrasts, direct between-sex differences, and pooled absolute rates with gamma intervals. |
| `scripts/period_age_sensitivity.R` | Contextual 2000–2020 analysis using adult age bands 20–59 and 60+, with paired contrasts for both sexes combined and separately. |
| `scripts/figures.R` | Optionally generate the three main incidence/composition plots and two supplementary diagnostic plots. |

The primary population is all malignant C23.9 tumors (D1). The SEER Rare Cancer Classification 15.1 / recode 60 category supplies a supporting analysis. Adenocarcinoma NOS (8140/3; D4) is a sensitivity population, not all adenocarcinomas. Numeric ranges 8010–8576 (D2) and 8140–8389 (D3) are used for coding audits, not clinically pure adenocarcinoma definitions. The populations overlap.

## Methods

Rates use the 2000 US age standard, normalized separately within ages 20–49 and 50+. Joint standardization fixes within-five-year-age-band ethnic population shares to 2000. The 2020 and pooled 2000–2023 person-year references provide sensitivity comparisons. Unknown-race cases are audited and excluded from paired rates because their population denominator is unavailable. Both comparators use the same known-ethnicity cases.

APCs are obtained from inverse-variance-weighted log-rate regressions; group-specific confidence intervals use the residual variance of a two-age-group interaction model. The contrast Δ is joint minus age-only APC, in percentage points/year, computed before rounding. It is not the annual change in a younger/older rate ratio.

Each calendar-window count generator is fitted only to its included years, with spline knots recalculated within that window. Positive fitted means allow nonzero simulated counts for observed zero cells. Both rates, variances and WLS weights are recalculated from the same count draw. Smooth-generator basic intervals are centered with Δfit, obtained by applying the entire estimator to fitted means. Plug-in Poisson percentile intervals and Jeffreys Gamma posterior credible intervals are distinct uncertainty analyses.

Annual-trend simulations use 10000 draws and seed 20261025; period decomposition uses a separate seed 20261026. Negative binomial NB2 sensitivity uses seed 261026, 10000 draws, and fixed fitted theta; theta uncertainty is not integrated. Sex-stratified analysis uses seed 2610261, 10000 paired draws, independent sex/window generators, and the common within-age 2000 ethnic shares pooled over both sexes. Own-sex references are secondary checks. Between-sex intervals use the replicate-wise difference between female and male Δ, not separate P values.

Period decomposition separates the absolute change into composition and within-stratum components using symmetric weights. Periods are 2000–2004 versus 2019–2023, 2000–2004 versus 2015–2019, and 2000 versus 2023. Contributions are descriptive, not causal or preventable fractions. Denominators and reference weights are fixed. No registry-level resampling, reporting-delay adjustment, mortality modeling or screening inference is performed. Subgroup comparisons are exploratory without multiplicity correction. Conditional independence of sex-specific draws does not capture shared registry dependence.

Generated outputs can include sparse cells. Keep them private until permitted sharing is confirmed. The in-memory tests check mathematical properties, not independent replication. Exact study-data reproduction requires the matching authorized exports. Consult the Zenodo archive metadata for its license. No separate software license file is included.

Pooled sex-specific absolute rates sum counts and population person-years within each age–ethnicity cell for 2000–2004 and 2019–2023 before standardization. Rates in `sex_pooled_rates.csv` are per 100000; multiply rates and interval limits by ten for presentation per million person-years. Joint rates use the common two-sex 2000 within-age ethnic shares. Gamma intervals use the largest coefficient over all eligible cells, including zero-count cells. No expected or preventable case counts are calculated. The input exports contain no registry dimension, so this workflow does not perform registry-exclusion or geography-standardization analyses.

Main tables 1–2 and supplementary tables S1–S7 and S9–S11 use the corresponding result CSVs; the coding audit and S8 ethnic count-model trends are generated under `private_inputs`. Main table 3 and S12–S13 use `sex_paired_trends.csv` and `sex_between_differences.csv`; pooled absolute rates in main table 3 and S14 use `sex_pooled_rates.csv`. S15 is a literature comparison and the STROBE file is a reporting checklist, not statistical outputs. S16 uses `period_age_sensitivity.csv`; its generator diagnostics and annual rates are also exported. The optional figure script generates the three main figures and supplementary Figures S1–S2. Document layout and reference formatting are not generated by this repository.

The contextual period/age sensitivity uses 2000–2020 and ages 20–59 versus 60+ under the current D1 definition and database. It uses a common two-sex 2000 ethnic reference, window-specific positive-mean generators, 10000 paired draws and seed 261031. It is not a parameter-matched replication of the 2024 HPB research letter: that letter used a November 2022 submission and an under-60 category, while its exact morphology code list, registry coverage, standard population and complete rate-session settings were not reported. The current inputs exclude ages below 20. These differences must remain explicit when interpreting the contextual results.

## Method sources

- [SEER Rare Cancer Classification](https://seer.cancer.gov/seerstat/variables/seer/raresiterecode/), category 15.1.
- [2000 US standard population](https://seer.cancer.gov/stdpopulations/stdpop.20ages.html).
- Symmetric decomposition: [Kitagawa (1955)](https://doi.org/10.1080/01621459.1955.10501299).
- Gamma intervals: [Fay and Feuer (1997)](https://doi.org/10.1002/(SICI)1097-0258(19970415)16:7%3C791::AID-SIM500%3E3.0.CO;2-%23).
