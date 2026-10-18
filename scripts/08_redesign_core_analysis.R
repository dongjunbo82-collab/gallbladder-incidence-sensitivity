options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(splines)
})

args <- commandArgs(TRUE)
root <- normalizePath(args[1], winslash="/")
raw_dir <- normalizePath(args[2], winslash="/")
table_dir <- file.path(root, "Analysis_Code_and_Source_Data", "derived_results")
qa_dir <- file.path(root, "Analysis_Code_and_Source_Data", "analysis_checks")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

read_export <- function(path) {
  read.csv(path, check.names = FALSE, na.strings = c("NA", ""),
           fileEncoding = "UTF-8-BOM")
}

num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))

reshape_p1 <- function(path, registry) {
  x <- read_export(path)
  year <- num(x[[1]])
  pop_cols <- grep("_Pop$", names(x), value = TRUE)
  count_cols <- grep("_Count$", names(x), value = TRUE)
  age_code <- sub("_Pop$", "", pop_cols)
  age_start <- num(sub("_.*$", "", age_code))
  age_end <- ifelse(age_start == 90, NA, age_start + 4)
  out <- do.call(rbind, lapply(seq_along(pop_cols), function(i) {
    data.frame(
      registry = registry,
      year = year,
      age_start = age_start[i],
      age_end = age_end[i],
      age_mid = ifelse(age_start[i] == 90, NA, (age_start[i] + age_end[i]) / 2),
      p1_count = num(x[[count_cols[i]]]),
      population = num(x[[pop_cols[i]]])
    )
  }))
  out[order(out$year, out$age_start), ]
}

parse_single_age <- function(x) {
  out <- num(sub(" years.*$", "", trimws(as.character(x))))
  out[grepl("^90\\+", trimws(as.character(x)))] <- 90
  out
}

age5_start <- function(age) {
  ifelse(age >= 90, 90, floor(age / 5) * 5)
}

definition_from_histology <- function(hist, definition) {
  if (definition == "D1_all_malignant") return(rep(TRUE, length(hist)))
  if (definition == "D2_carcinoma_8010_8576") return(hist >= 8010 & hist <= 8576)
  if (definition == "D3_adenocarcinoma_8140_8389") return(hist >= 8140 & hist <= 8389)
  if (definition == "D4_strict_8140_3") return(hist == 8140)
  stop("Unknown definition")
}

definitions <- c(
  "D1_all_malignant",
  "D2_carcinoma_8010_8576",
  "D3_adenocarcinoma_8140_8389",
  "D4_strict_8140_3"
)

# Population denominators and case-level counts ---------------------------
p1_17 <- reshape_p1(
  file.path(raw_dir, "P1_GB_AllMalignant_YearAge5_Rates_SEER17_2000-2023.csv"),
  "SEER17"
)
p1_8 <- reshape_p1(
  file.path(raw_dir, "V8_P1_GB_AllMalignant_YearAge5_Rates_SEER8_2000-2023.csv"),
  "SEER8"
)

c2 <- read_export(file.path(raw_dir, "C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt"))
c2$year <- num(c2[["Year of diagnosis"]])
c2$age <- parse_single_age(c2[["Age recode with single ages and 90+"]])
c2$age_start <- age5_start(c2$age)
c2$histology <- num(c2[["Histologic Type ICD-O-3"]])
c2 <- c2[is.finite(c2$year) & is.finite(c2$age) & is.finite(c2$histology) &
           c2$year >= 2000 & c2$year <= 2023 & c2$age >= 20, ]

count_list <- lapply(definitions, function(definition) {
  z <- c2[definition_from_histology(c2$histology, definition), ]
  a <- aggregate(rep(1, nrow(z)), by = list(year = z$year, age_start = z$age_start), sum)
  names(a)[3] <- "count"
  a$definition <- definition
  a
})
case_counts <- do.call(rbind, count_list)

grid17 <- do.call(rbind, lapply(definitions, function(definition) {
  z <- p1_17
  z$definition <- definition
  z
}))
incidence_cells <- merge(grid17, case_counts,
                         by = c("year", "age_start", "definition"), all.x = TRUE)
incidence_cells$count[is.na(incidence_cells$count)] <- 0
incidence_cells <- incidence_cells[order(incidence_cells$definition,
                                         incidence_cells$year,
                                         incidence_cells$age_start), ]

# Core data QA ------------------------------------------------------------
d1_cells <- incidence_cells[incidence_cells$definition == "D1_all_malignant", ]
qa <- data.frame(
  check = c(
    "C2 D1 total equals P1 total",
    "D1 counts match P1 in every year-age cell",
    "All four definitions have 360 year-age cells",
    "No nonpositive population denominators",
    "Definitions are nested in every year-age cell",
    "SEER8 P1 has 360 year-age cells"
  ),
  observed = c(
    paste(sum(d1_cells$count), sum(p1_17$p1_count), sep = "/"),
    sum(d1_cells$count != d1_cells$p1_count),
    paste(as.integer(table(incidence_cells$definition)), collapse = "/"),
    sum(incidence_cells$population <= 0 | !is.finite(incidence_cells$population)),
    NA_character_,
    nrow(p1_8)
  ),
  expected = c("24509/24509", "0", "360/360/360/360", "0", "0", "360"),
  pass = FALSE
)

wide_counts <- reshape(incidence_cells[, c("year", "age_start", "definition", "count")],
                       idvar = c("year", "age_start"), timevar = "definition",
                       direction = "wide")
nested_viol <- sum(
  wide_counts$count.D1_all_malignant < wide_counts$count.D2_carcinoma_8010_8576 |
    wide_counts$count.D2_carcinoma_8010_8576 < wide_counts$count.D3_adenocarcinoma_8140_8389 |
    wide_counts$count.D3_adenocarcinoma_8140_8389 < wide_counts$count.D4_strict_8140_3
)
qa$observed[5] <- nested_viol
qa$pass <- as.character(qa$observed) == qa$expected

# Continuous-age incidence models ----------------------------------------
contrast_from_fit <- function(fit, new1, new2) {
  tt <- delete.response(terms(fit))
  x1 <- model.matrix(tt, new1, contrasts.arg = fit$contrasts)
  x2 <- model.matrix(tt, new2, contrasts.arg = fit$contrasts)
  x1 <- x1[, names(coef(fit)), drop = FALSE]
  x2 <- x2[, names(coef(fit)), drop = FALSE]
  d <- as.numeric(x2 - x1)
  b <- sum(d * coef(fit))
  se <- sqrt(as.numeric(t(d) %*% vcov(fit) %*% d))
  c(estimate = b, se = se)
}

fit_continuous_age <- function(z, definition) {
  # The 90+ denominator is an open interval and is therefore retained in
  # grouped-rate summaries but excluded from the continuous-age spline.
  z <- z[z$age_start < 90 & is.finite(z$age_mid), ]
  z$year_c <- z$year - 2011.5
  full <- glm(count ~ ns(age_mid, df = 4) * year_c + offset(log(population)),
              family = quasipoisson(), data = z)
  reduced <- glm(count ~ ns(age_mid, df = 4) + year_c + offset(log(population)),
                 family = quasipoisson(), data = z)
  interaction_test <- anova(reduced, full, test = "F")
  ages <- c(30, 40, 50, 60, 70, 80)
  slopes <- do.call(rbind, lapply(ages, function(a) {
    n1 <- data.frame(age_mid = a, year_c = 0, population = 100000)
    n2 <- data.frame(age_mid = a, year_c = 1, population = 100000)
    cs <- contrast_from_fit(full, n1, n2)
    data.frame(
      definition = definition,
      age = a,
      apc = 100 * (exp(cs["estimate"]) - 1),
      lower_95 = 100 * (exp(cs["estimate"] - 1.96 * cs["se"]) - 1),
      upper_95 = 100 * (exp(cs["estimate"] + 1.96 * cs["se"]) - 1),
      p_value = 2 * pnorm(abs(cs["estimate"] / cs["se"]), lower.tail = FALSE),
      age_by_year_global_p = interaction_test$`Pr(>F)`[2],
      dispersion = summary(full)$dispersion
    )
  }))
  years <- c(2000, 2010, 2023)
  curve_ages <- seq(22, 87, 5)
  curves <- do.call(rbind, lapply(years, function(y) {
    nd <- data.frame(age_mid = curve_ages, year_c = y - 2011.5, population = 100000)
    pr <- predict(full, newdata = nd, type = "link", se.fit = TRUE)
    data.frame(
      definition = definition,
      year = y,
      age_mid = curve_ages,
      predicted_rate = exp(pr$fit),
      lower_95 = exp(pr$fit - 1.96 * pr$se.fit),
      upper_95 = exp(pr$fit + 1.96 * pr$se.fit)
    )
  }))
  list(fit = full, slopes = slopes, curves = curves)
}

continuous_fits <- lapply(definitions, function(definition) {
  fit_continuous_age(incidence_cells[incidence_cells$definition == definition, ], definition)
})
continuous_slopes <- do.call(rbind, lapply(continuous_fits, `[[`, "slopes"))
continuous_curves <- do.call(rbind, lapply(continuous_fits, `[[`, "curves"))

# Annual official-rate summaries for each registry and definition ---------
rate_files_17 <- c(
  D1_all_malignant = "R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER17_2000-2023.txt",
  D2_carcinoma_8010_8576 = "R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER17_2000-2023.txt",
  D3_adenocarcinoma_8140_8389 = "R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER17_2000-2023.txt",
  D4_strict_8140_3 = "R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER17_2000-2023.txt"
)
rate_files_8 <- c(
  D1_all_malignant = "V8_R6_GB_D1_AllMalignant_AgeGroup_AAR_SEER8_2000-2023.txt",
  D2_carcinoma_8010_8576 = "V8_R7_GB_D2_Carcinoma8010-8576_AgeGroup_AAR_SEER8_2000-2023.txt",
  D3_adenocarcinoma_8140_8389 = "V8_R8_GB_D3_Adenocarcinoma8140-8389_AgeGroup_AAR_SEER8_2000-2023.txt",
  D4_strict_8140_3 = "V8_R9_GB_D4_Strict8140-3_AgeGroup_AAR_SEER8_2000-2023.txt"
)

read_age_group_rates <- function(path, registry, definition) {
  x <- read_export(path)
  do.call(rbind, lapply(c("20_49", "50"), function(prefix) {
    data.frame(
      registry = registry,
      definition = definition,
      group = ifelse(prefix == "20_49", "20-49", "50+"),
      year = num(x[[1]]),
      rate = num(x[[paste0(prefix, "_Rate")]]),
      se = num(x[[paste0(prefix, "_SE")]]),
      count = num(x[[paste0(prefix, "_Count")]]),
      population = num(x[[paste0(prefix, "_Pop")]])
    )
  }))
}

all_rates <- do.call(rbind, c(
  Map(function(f, d) read_age_group_rates(file.path(raw_dir, f), "SEER17", d),
      rate_files_17, names(rate_files_17)),
  Map(function(f, d) read_age_group_rates(file.path(raw_dir, f), "SEER8", d),
      rate_files_8, names(rate_files_8))
))

wls_apc <- function(z) {
  z <- z[is.finite(z$rate) & is.finite(z$se) & z$rate > 0 & z$se > 0, ]
  fit <- lm(log(rate) ~ year, weights = 1 / (se / rate)^2, data = z)
  b <- coef(fit)["year"]
  s <- sqrt(vcov(fit)["year", "year"])
  data.frame(
    apc = 100 * (exp(b) - 1),
    lower_95 = 100 * (exp(b - qt(.975, df.residual(fit)) * s) - 1),
    upper_95 = 100 * (exp(b + qt(.975, df.residual(fit)) * s) - 1),
    p_value = 2 * pt(abs(b / s), df.residual(fit), lower.tail = FALSE)
  )
}

trend_summary <- do.call(rbind, lapply(split(all_rates,
  interaction(all_rates$registry, all_rates$definition, all_rates$group, drop = TRUE)), function(z) {
  cbind(z[1, c("registry", "definition", "group")], wls_apc(z))
}))
rownames(trend_summary) <- NULL

trend_interactions <- do.call(rbind, lapply(split(all_rates,
  interaction(all_rates$registry, all_rates$definition, drop = TRUE)), function(z) {
  z$group <- relevel(factor(z$group), ref = "50+")
  fit <- lm(log(rate) ~ year * group, weights = 1 / (se / rate)^2, data = z)
  term <- grep("year:group", names(coef(fit)), value = TRUE)[1]
  b <- coef(fit)[term]
  s <- sqrt(vcov(fit)[term, term])
  data.frame(
    registry = z$registry[1], definition = z$definition[1],
    differential_apc = 100 * (exp(b) - 1),
    lower_95 = 100 * (exp(b - qt(.975, df.residual(fit)) * s) - 1),
    upper_95 = 100 * (exp(b + qt(.975, df.residual(fit)) * s) - 1),
    p_interaction = 2 * pt(abs(b / s), df.residual(fit), lower.tail = FALSE)
  )
}))
rownames(trend_interactions) <- NULL

# Standardized cutoff and calendar-window analyses are implemented in script 11.

# Case-definition composition mechanisms ---------------------------------
c2$age_group <- ifelse(c2$age < 50, "20-49", "50+")
c2$period <- cut(c2$year, breaks = c(1999, 2004, 2009, 2014, 2019, 2023),
                 labels = c("2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2023"))
c2$definition_category <- ifelse(c2$histology == 8140, "D4 strict 8140/3",
  ifelse(c2$histology >= 8140 & c2$histology <= 8389, "D3 other adenocarcinoma",
    ifelse(c2$histology >= 8010 & c2$histology <= 8576, "D2 other carcinoma",
           "Other malignant histology")))
composition_counts <- as.data.frame(xtabs(~ age_group + period + definition_category, data = c2))
composition_counts$percent_within_age_period <- ave(composition_counts$Freq,
  interaction(composition_counts$age_group, composition_counts$period),
  FUN = function(x) 100 * x / sum(x))

c2$strict <- as.integer(c2$histology == 8140)
c2$year_c <- c2$year - 2011.5
c2$young <- factor(ifelse(c2$age < 50, "20-49", "50+"), levels = c("50+", "20-49"))
c2$microscopic <- grepl("histology|cytology|microscopic|immunophenotyping|genetic",
                        c2[["Diagnostic Confirmation"]], ignore.case = TRUE)
c2$hospital_report <- grepl("Hospital|surgery center|oncology center",
                            c2[["Type of Reporting Source"]], ignore.case = TRUE)
composition_unadjusted <- glm(strict ~ year_c * young + ns(age, df = 4),
                              family = binomial(), data = c2)
composition_adjusted <- glm(strict ~ year_c * young + ns(age, df = 4) +
                              microscopic + hospital_report,
                            family = binomial(), data = c2)

extract_interaction <- function(fit, model) {
  term <- grep("year_c:young", names(coef(fit)), value = TRUE)[1]
  b <- coef(fit)[term]
  s <- sqrt(vcov(fit)[term, term])
  data.frame(
    model = model,
    odds_ratio_per_year = exp(b),
    lower_95 = exp(b - 1.96 * s),
    upper_95 = exp(b + 1.96 * s),
    p_interaction = 2 * pnorm(abs(b / s), lower.tail = FALSE),
    n = nobs(fit)
  )
}
composition_models <- rbind(
  extract_interaction(composition_unadjusted, "Age-year model"),
  extract_interaction(composition_adjusted, "Adjusted for confirmation and reporting source")
)

# Write outputs -----------------------------------------------------------
write.csv(qa, file.path(qa_dir, "redesign_core_data_QA_20261011.csv"), row.names = FALSE)
write.csv(incidence_cells, file.path(table_dir, "incidence_year_age_definition_cells.csv"), row.names = FALSE)
write.csv(continuous_slopes, file.path(table_dir, "continuous_age_apc_estimates.csv"), row.names = FALSE)
write.csv(continuous_curves, file.path(table_dir, "continuous_age_predicted_rates.csv"), row.names = FALSE)
write.csv(all_rates, file.path(table_dir, "official_age_group_rates_seer17_seer8.csv"), row.names = FALSE)
write.csv(trend_summary, file.path(table_dir, "loglinear_apc_by_definition_registry.csv"), row.names = FALSE)
write.csv(trend_interactions, file.path(table_dir, "age_group_trend_interactions_by_definition_registry.csv"), row.names = FALSE)
# Script 11 overwrites all SEER17 grouped summaries using a single standardized
# rate framework, and is the authoritative source for cutoff sensitivity.
write.csv(composition_counts, file.path(table_dir, "definition_composition_by_age_period.csv"), row.names = FALSE)
write.csv(composition_models, file.path(table_dir, "definition_composition_interaction_models.csv"), row.names = FALSE)

audit <- c(
  "GBAC redesign core incidence analysis",
  paste("Run time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Core data QA:", ifelse(all(qa$pass), "PASS", "FAIL")),
  paste("SEER17 cases by definition:", paste(tapply(incidence_cells$count,
    incidence_cells$definition, sum), collapse = ", ")),
  paste("Continuous-age global interaction p values:",
        paste(sprintf("%s=%.4g", continuous_slopes$definition[!duplicated(continuous_slopes$definition)],
                      continuous_slopes$age_by_year_global_p[!duplicated(continuous_slopes$definition)]),
              collapse = "; ")),
  "Trend estimates from R6-R9 are log-linear WLS summaries, not official NCI Joinpoint estimates.",
  "SEER8 is an overlapping-registry sensitivity analysis, not independent validation.",
  "No manuscript text was modified by this analysis."
)
writeLines(audit, file.path(qa_dir, "redesign_core_analysis_audit_20261011.txt"), useBytes = TRUE)

cat(paste(audit, collapse = "\n"), "\n")
