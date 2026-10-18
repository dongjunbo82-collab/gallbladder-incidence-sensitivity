# Run the complete analysis without including any study data in this repository.
args <- commandArgs(TRUE)
if (length(args) < 2L || any(args %in% c('-h', '--help'))) {
  cat('Usage: Rscript run_analysis.R OUTPUT_DIRECTORY RAW_DIRECTORY [--no-figures]\n')
  quit(status = if (any(args %in% c('-h', '--help'))) 0L else 1L)
}
file_arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop('Run this entry point with Rscript.')
project <- dirname(normalizePath(sub('^--file=', '', file_arg), winslash = '/'))
raw <- normalizePath(args[2], winslash = '/', mustWork = TRUE)
out_parent <- normalizePath(dirname(path.expand(args[1])), winslash = '/', mustWork = TRUE)
out <- file.path(out_parent, basename(args[1]))
if (out %in% c(raw, project) || startsWith(out, paste0(raw, '/')) ||
    startsWith(out, paste0(project, '/'))) stop('Output must be outside the raw and code directories.')
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) {
  stop('Choose a new or empty output directory; existing outputs are not overwritten.')
}
required_packages <- c('splines', 'survival', 'cmprsk')
figures <- !('--no-figures' %in% args)
if (figures) required_packages <- c(required_packages, 'ggplot2', 'patchwork', 'scales', 'ragg')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) stop('Required packages unavailable: ', paste(missing_packages, collapse = ', '), '. No installation was attempted.')
required_files <- c(
  'C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt',
  'P1_GB_AllMalignant_YearAge5_Rates_SEER17_2000-2023.csv',
  'V8_P1_GB_AllMalignant_YearAge5_Rates_SEER8_2000-2023.csv',
  'P2_GB_AllMalignant_YearAge5_RaceEthnicity_SEER17_2000-2023.txt',
  'R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER17_2000-2023.txt',
  'R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER17_2000-2023.txt',
  'R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER17_2000-2023.txt',
  'R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER17_2000-2023.txt',
  'V8_R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER8_2000-2023.txt',
  'V8_R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER8_2000-2023.txt',
  'V8_R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER8_2000-2023.txt',
  'V8_R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER8_2000-2023.txt',
  'R2_GBAC8140_EO20-49_AnnualRates_bySex_SEER17_2000-2023.csv',
  'R3_GBAC8140_EO20-49_5YRates_byRaceEthnicity_SEER17_2000-2023.csv',
  'R4b_GBAC8140_EO20-49_StageSpecificRates_SEER17_2004-2023.csv'
)
missing_files <- required_files[!file.exists(file.path(raw, required_files))]
if (length(missing_files)) stop('Missing authorized local exports:\n', paste(missing_files, collapse = '\n'))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script_out <- file.path(out, 'Analysis_Code_and_Source_Data', 'scripts')
dir.create(script_out, recursive = TRUE, showWarnings = FALSE)
scripts <- list.files(file.path(project, 'scripts'), pattern = '[.]R$', full.names = TRUE)
stopifnot(all(file.copy(scripts, script_out)))
log_dir <- file.path(out, 'execution_logs'); dir.create(log_dir)
writeLines(capture.output(sessionInfo()), file.path(log_dir, 'environment.txt'))
order <- c('08_redesign_core_analysis.R', '09_competing_risk_analysis.R',
           '11_revision_analysis.R', '12_time_varying_mortality.R',
           '13_ethnicity_analysis.R', '14_reframe_analysis.R', '15_robustness_analysis.R')
if (figures) order <- c(order, '10_manuscript_figures.R')
rscript <- file.path(R.home('bin'), if (.Platform$OS.type == 'windows') 'Rscript.exe' else 'Rscript')
for (name in order) {
  cat('Running ', name, '\n', sep = ''); flush.console()
  second <- if (name == '10_manuscript_figures.R') file.path(out, 'internal_previews') else raw
  log <- file.path(log_dir, paste0(name, '.log'))
  status <- system2(rscript, shQuote(c(file.path(script_out, name), out, second)), stdout = log, stderr = log)
  if (status != 0L) stop('Analysis failed in ', name, '; inspect its local execution log.')
}
cat('Completed. Generated outputs remain local, outside this code repository.\n')
