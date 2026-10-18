# Local input requirements

This document describes filenames and field names only. It contains no records,
sample rows, extracted rates or result tables. Keep all input files outside the
repository and never commit generated outputs.

## Input filenames

The runner requires 15 files from the same November 2025 SEER submission:

```text
C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt
P1_GB_AllMalignant_YearAge5_Rates_SEER17_2000-2023.csv
V8_P1_GB_AllMalignant_YearAge5_Rates_SEER8_2000-2023.csv
P2_GB_AllMalignant_YearAge5_RaceEthnicity_SEER17_2000-2023.txt
R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER17_2000-2023.txt
R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER17_2000-2023.txt
R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER17_2000-2023.txt
R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER17_2000-2023.txt
V8_R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER8_2000-2023.txt
V8_R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER8_2000-2023.txt
V8_R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER8_2000-2023.txt
V8_R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER8_2000-2023.txt
R2_GBAC8140_EO20-49_AnnualRates_bySex_SEER17_2000-2023.csv
R3_GBAC8140_EO20-49_5YRates_byRaceEthnicity_SEER17_2000-2023.csv
R4b_GBAC8140_EO20-49_StageSpecificRates_SEER17_2004-2023.csv
```

These are comma-delimited exports with preserved column names, despite some
`.txt` extensions. Retain exact labels, zero cells and the matching dictionaries
locally. Do not replace unavailable or suppressed values with zero. The fixed
shape/count assertions in the scripts are validation constants, not distributed
datasets. Exports with different sessions or database releases require review of
those assertions and variable definitions before use.

## Case-listing fields

C2 covers malignant primary site C23.9, diagnosis years 2000-2023, ages 20+.
Required labeled columns include:

- Year of diagnosis
- Age recode with single ages and 90+
- Histologic Type ICD-O-3
- ICD-O-3 Hist/behav, malignant
- Sex
- Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)
- Diagnostic Confirmation
- Type of Reporting Source
- Survival months
- SEER cause-specific death classification
- SEER other cause of death classification
- Combined Summary Stage with Expanded Regional Codes (2004+)
- Grade Recode (thru 2017)
- Grade Pathological (2018+)
- Grade Clinical (2018+)

The analysis does not need personal identifiers. Do not add names, addresses or
other identifying fields when obtaining or preparing the local exports.

## Population and rate fields

P1 files use an annual year column and wide age-band columns ending in `_Count`
and `_Pop`. Age-band column prefixes begin with the lower age bound; the bands
cover ages 20-24 through 85-89 and 90+. R6-R9 and V8 counterparts have an annual
year column and `20_49_Rate`, `20_49_SE`, `20_49_Count`, `20_49_Pop`, `50_Rate`,
`50_SE`, `50_Count`, `50_Pop` columns.

P2 is long-form and includes `Year_Annual_2000_2023`,
`Age recode with <1 year olds and 90+`, the five-category race/origin field,
`Count` and `Population`. Unknown race has no population denominator and is
audited separately. Structural excluded-age placeholders must have zero count
and population. Known-ethnicity population cells must be positive.

R2 includes `Sex`, `Year of diagnosis`, `Age-Adjusted Rate/Trend`,
`Standard Error/P-Value`, native trend confidence limits and APC rows.
R4b includes `STAGE_LRD`, `PERIOD_STAGE_5Y`, `Age-Adjusted Rate`,
`Standard Error`, `Lower Confidence Interval`, `Upper Confidence Interval`,
`Count` and `Population`. R3 is retained as a descriptive input summary.

Stage labels are Localized, Regional, Distant and Unknown. The four stage-rate
periods start in 2004. Check native dictionaries before substituting custom recodes.
