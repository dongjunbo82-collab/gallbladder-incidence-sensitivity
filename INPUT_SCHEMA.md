# Authorized input exports

Supply three labeled CSV-format text exports in a private input directory. Despite the `.txt` extension, the verified exports are comma-delimited with a header row and are read with `read.csv(check.names=FALSE)`. Remove thousands separators in exported numeric fields. Retain zero-count cells and do not suppress cells in analytical exports. No input or result data are distributed.

## Shared eligibility

SEER Research Data, 17 Registries, November 2025 submission, diagnoses 2000–2023. Primary site C23.9 (gallbladder), malignant behavior /3, ages ≥20 including 90+, both sexes and all tumor sequences. No histology restriction in source sessions. Rate sessions must use the matching SEER population geography. These rates are not delay-adjusted.

The case listing is assumed to have applied these eligibility restrictions. If available fields include site or behavior, check these restrictions before running. The workflow cannot verify an export's database identity from its cells alone. Keep the original session and dictionary files locally for provenance.

## C2 case listing

Filename: `C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt`

Required labeled fields:

- `Age recode with single ages and 90+`
- `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`
- `Year of diagnosis`
- `Histologic Type ICD-O-3`

One row represents one eligible tumor. No identifiers, outcome or treatment columns are needed. The analysis uses the four named fields only.

## P2 rate cells

Filename: `P2_GB_AllMalignant_YearAge5_RaceEthnicity_SEER17_2000-2023.txt`

Required fields: `Year_Annual_2000_2023`, `Age recode with <1 year olds and 90+`, `Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)`, `Count`, `Population`.

Use annual years, five-year age bands 20–24 through 85–89 and 90+, and the five known race/origin categories plus `Non-Hispanic Unknown Race`. Known categories have positive denominators; unknown race has zero population. Display placeholders below age 20 may remain only if count and population are both zero. Eligible P2 cells number 2160, of which 1800 have known race. C2 counts must exactly equal P2 counts in every eligible cell.

## P3 sex-specific rate cells

Filename: `P3_GB_AllMalignant_YearAge5_Sex_RaceEthnicity_SEER17_2000-2023.txt`

Required fields are those of P2 plus `Sex`, labeled `Male` and `Female`. Total-sex display rows may be present and are excluded. Male plus female counts and populations must exactly equal P2 in each race/year/age cell. There are 3600 known-race male/female eligible cells. Unknown race is audited separately. The `Rate` field is not used; rates are recalculated from counts and denominators.

## Generated private files

`prepare_inputs.R` creates all denominator cells, classification lists, reference shares, WLS summaries and pooled decomposition cells required by the remaining scripts. These files are computational intermediates, not public inputs. Results and logs are written under the output directory. Do not upload them with the code.
