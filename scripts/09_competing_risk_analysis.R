options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
  library(splines)
})

args <- commandArgs(TRUE)
root <- normalizePath(args[1], winslash="/")
raw_file <- file.path(normalizePath(args[2], winslash="/"), "C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt")
table_dir <- file.path(root, "Analysis_Code_and_Source_Data", "derived_results")
qa_dir <- file.path(root, "Analysis_Code_and_Source_Data", "analysis_checks")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(raw_file, check.names = FALSE, na.strings = c("NA", ""),
              fileEncoding = "UTF-8-BOM")
num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))

d$year <- num(d[["Year of diagnosis"]])
d$age <- num(sub(" years.*$", "", trimws(d[["Age recode with single ages and 90+"]])))
d$age[grepl("^90\\+", trimws(d[["Age recode with single ages and 90+"]]))] <- 90
d$histology <- num(d[["Histologic Type ICD-O-3"]])
d$time_months <- num(d[["Survival months"]])
d$age_group <- factor(ifelse(d$age < 50, "20-49", "50+"), levels = c("50+", "20-49"))
d$sex <- factor(d$Sex)
d$race <- factor(d[["Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)"]])
d$report_source <- d[["Type of Reporting Source"]]
d$cancer_death_class <- d[["SEER cause-specific death classification"]]
d$other_death_class <- d[["SEER other cause of death classification"]]

d$stage <- ifelse(grepl("^Localized", d[["Combined Summary Stage with Expanded Regional Codes (2004+)"]]), "Localized",
  ifelse(grepl("^Regional", d[["Combined Summary Stage with Expanded Regional Codes (2004+)"]]), "Regional",
    ifelse(grepl("^Distant", d[["Combined Summary Stage with Expanded Regional Codes (2004+)"]]), "Distant", "Unknown")))
d$stage <- factor(d$stage, levels = c("Localized", "Regional", "Distant", "Unknown"))

grade_old <- d[["Grade Recode (thru 2017)"]]
grade_path <- d[["Grade Pathological (2018+)"]]
grade_clin <- d[["Grade Clinical (2018+)"]]
grade <- rep("Unknown", nrow(d))
grade[grepl("Grade I$|Well differentiated", grade_old)] <- "I"
grade[grepl("Grade II$|Moderately differentiated", grade_old)] <- "II"
grade[grepl("Grade III|Grade IV|Poorly differentiated|Undifferentiated", grade_old)] <- "III-IV"
g2018 <- ifelse(grade_path %in% c("1", "2", "3"), grade_path,
                ifelse(grade_clin %in% c("1", "2", "3"), grade_clin, NA))
grade[d$year >= 2018 & g2018 == "1"] <- "I"
grade[d$year >= 2018 & g2018 == "2"] <- "II"
grade[d$year >= 2018 & g2018 == "3"] <- "III-IV"
d$grade <- factor(grade, levels = c("I", "II", "III-IV", "Unknown"))

d$status <- ifelse(d$cancer_death_class == "Dead (attributable to this cancer dx)", 1,
  ifelse(d$other_death_class == "Dead (attributable to causes other than this cancer dx)", 2, 0))
d$unknown_cod <- grepl("missing/unknown COD", d$cancer_death_class, fixed = TRUE)
d$dco_autopsy <- grepl("Death certificate only|Autopsy only", d$report_source)

step_estimate <- function(time, estimate, variance, target) {
  idx <- max(which(time <= target), 0)
  if (idx == 0) return(c(est = 0, se = 0))
  c(est = estimate[idx], se = sqrt(variance[idx]))
}

crude_cif <- function(z, definition, window, horizon) {
  z$time_h <- pmin(z$time_months, horizon)
  z$status_h <- ifelse(z$time_months > horizon, 0, z$status)
  ci <- cuminc(z$time_h, z$status_h, group = z$age_group, cencode = 0)
  targets <- sort(unique(c(12, 36, horizon)))
  out <- do.call(rbind, lapply(levels(z$age_group), function(g) {
    obj <- ci[[paste(g, "1")]]
    do.call(rbind, lapply(targets, function(t) {
      e <- step_estimate(obj$time, obj$est, obj$var, t)
      data.frame(
        definition = definition, window = window, age_group = g,
        months = t, cif = e["est"],
        lower_95 = max(0, e["est"] - 1.96 * e["se"]),
        upper_95 = min(1, e["est"] + 1.96 * e["se"]),
        n = sum(z$age_group == g),
        cancer_deaths = sum(z$age_group == g & z$status_h == 1),
        other_deaths = sum(z$age_group == g & z$status_h == 2)
      )
    }))
  }))
  gray_p <- ci$Tests["1", "pv"]
  out$gray_p <- gray_p
  out
}

model_results <- function(z, definition, window, horizon) {
  z$time_h <- pmin(z$time_months, horizon)
  z$status_h <- ifelse(z$time_months > horizon, 0, z$status)
  z$year_c <- z$year - mean(z$year)
  mm <- model.matrix(~ age_group + sex + race + year_c + stage + grade, data = z)[, -1, drop = FALSE]
  fg <- crr(z$time_h, z$status_h, cov1 = mm, failcode = 1, cencode = 0,
            variance = TRUE)
  fg_term <- grep("^age_group", names(fg$coef), value = TRUE)[1]
  fg_b <- fg$coef[fg_term]
  fg_se <- sqrt(diag(fg$var))[match(fg_term, names(fg$coef))]

  cs <- coxph(Surv(time_h, status_h == 1) ~ age_group + sex + race + year_c + stage + grade,
              data = z, ties = "efron")
  cs_term <- grep("^age_group", names(coef(cs)), value = TRUE)[1]
  cs_b <- coef(cs)[cs_term]
  cs_se <- sqrt(vcov(cs)[cs_term, cs_term])

  data.frame(
    definition = definition, window = window, horizon_months = horizon,
    comparison = "20-49 vs 50+",
    fine_gray_shr = exp(fg_b),
    fine_gray_lower_95 = exp(fg_b - 1.96 * fg_se),
    fine_gray_upper_95 = exp(fg_b + 1.96 * fg_se),
    fine_gray_p = 2 * pnorm(abs(fg_b / fg_se), lower.tail = FALSE),
    cause_specific_hr = exp(cs_b),
    cause_specific_lower_95 = exp(cs_b - 1.96 * cs_se),
    cause_specific_upper_95 = exp(cs_b + 1.96 * cs_se),
    cause_specific_p = 2 * pnorm(abs(cs_b / cs_se), lower.tail = FALSE),
    n = nrow(z), cancer_deaths = sum(z$status_h == 1),
    other_deaths = sum(z$status_h == 2)
  )
}

cohort_specs <- list(
  list(definition = "D4_strict_8140_3", hist = function(h) h == 8140,
       start = 2004, end = 2018, horizon = 60, window = "2004-2018_5year"),
  list(definition = "D4_strict_8140_3", hist = function(h) h == 8140,
       start = 2004, end = 2020, horizon = 36, window = "2004-2020_3year"),
  list(definition = "D3_adenocarcinoma_8140_8389", hist = function(h) h >= 8140 & h <= 8389,
       start = 2004, end = 2018, horizon = 60, window = "2004-2018_5year"),
  list(definition = "D3_adenocarcinoma_8140_8389", hist = function(h) h >= 8140 & h <= 8389,
       start = 2004, end = 2020, horizon = 36, window = "2004-2020_3year")
)

cohorts <- lapply(cohort_specs, function(s) {
  z <- d[s$hist(d$histology) & d$year >= s$start & d$year <= s$end &
           !d$dco_autopsy & !d$unknown_cod & is.finite(d$time_months), ]
  z$time_months[z$time_months == 0] <- 0.5
  list(spec = s, data = droplevels(z))
})

cohort_qa <- do.call(rbind, lapply(cohorts, function(x) {
  z <- x$data; s <- x$spec
  data.frame(
    definition = s$definition, window = s$window,
    n = nrow(z), young = sum(z$age_group == "20-49"), older = sum(z$age_group == "50+"),
    cancer_deaths = sum(z$status == 1), other_deaths = sum(z$status == 2),
    censored = sum(z$status == 0), zero_months_shifted = sum(z$time_months == 0.5),
    unknown_stage = sum(z$stage == "Unknown"), unknown_grade = sum(z$grade == "Unknown")
  )
}))

primary <- cohorts[[1]]$data
characteristic_rows <- function(var, label) {
  tab <- table(primary$age_group, primary[[var]], useNA = "ifany")
  lev <- colnames(tab)
  do.call(rbind, lapply(lev, function(lv) {
    data.frame(
      characteristic = label,
      level = lv,
      older_n = unname(tab["50+", lv]),
      older_pct = 100 * unname(tab["50+", lv]) / sum(primary$age_group == "50+"),
      young_n = unname(tab["20-49", lv]),
      young_pct = 100 * unname(tab["20-49", lv]) / sum(primary$age_group == "20-49")
    )
  }))
}

primary_characteristics <- rbind(
  characteristic_rows("sex", "Sex"),
  characteristic_rows("race", "Race and ethnicity"),
  characteristic_rows("stage", "Summary stage"),
  characteristic_rows("grade", "Grade")
)

cif_results <- do.call(rbind, lapply(cohorts, function(x) {
  s <- x$spec
  crude_cif(x$data, s$definition, s$window, s$horizon)
}))

adjusted_results <- do.call(rbind, lapply(cohorts, function(x) {
  s <- x$spec
  model_results(x$data, s$definition, s$window, s$horizon)
}))

write.csv(cohort_qa, file.path(qa_dir, "competing_risk_cohort_QA_20261011.csv"), row.names = FALSE)
write.csv(primary_characteristics, file.path(table_dir, "primary_survival_cohort_characteristics.csv"), row.names = FALSE)
write.csv(cif_results, file.path(table_dir, "crude_cancer_death_cif.csv"), row.names = FALSE)
write.csv(adjusted_results, file.path(table_dir, "adjusted_competing_risk_models.csv"), row.names = FALSE)

audit <- c(
  "GBAC competing-risk analysis",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "D4 denotes adenocarcinoma NOS; D3 denotes the 8140-8389 numeric morphology range.",
  "DCO/autopsy, unknown cause of death, and missing survival time were excluded.",
  "Zero-month survival was shifted to 0.5 month before time-to-event modeling.",
  "Crude cancer-death cumulative incidence accounts for other-cause death as a competing event.",
  "Adjusted associations include sex, race/ethnicity, diagnosis year, summary stage, and harmonized grade.",
  "Fine-Gray and cause-specific Cox estimates are associational, not causal treatment effects.",
  capture.output(print(adjusted_results, row.names = FALSE))
)
writeLines(audit, file.path(qa_dir, "competing_risk_analysis_audit_20261011.txt"), useBytes = TRUE)
cat(paste(audit, collapse = "\n"), "\n")
