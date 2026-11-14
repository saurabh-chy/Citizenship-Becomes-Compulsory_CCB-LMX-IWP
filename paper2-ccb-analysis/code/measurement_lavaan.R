# ============================================================================
# CCB TWO-PHASE STUDY: COMPLETE R/LAVAAN ANALYSIS
# ---------------------------------------------------------------------------
# This script produces:
#   1. Ordinal six-factor CFA versus one-factor CFA (WLSMV)
#   2. Standardised item loadings
#   3. Longitudinal ordinal LMX measurement invariance (configural, metric,
#      and scalar/threshold)
#   4. First-order SEM: Task, Contextual performance, and reverse-CWB as
#      separate latent performance outcomes
#   5. Second-order SEM: IWP as a higher-order factor represented by Task,
#      Contextual performance, and reverse-CWB
#   6. AMOS-style full-measurement and simplified structural SEM diagrams
#   7. Fit tables, standardised estimates, and bootstrap confidence intervals
#
# Design note:
# CCB T1 -> LMX T2 is the only baseline-adjusted temporal association.
# The performance paths are two-phase associations because LMX T2 and the
# performance outcomes were measured in Phase 2.
# ============================================================================
 
# -------------------------- 0. USER SETTINGS -------------------------------
# Usage: Rscript measurement_lavaan.R path/to/CCB_Two_Phase_Data.csv [output_folder]
args <- commandArgs(trailingOnly = TRUE)
DATA_PATH <- if (length(args) >= 1) args[1] else "CCB_Two_Phase_Data.csv"
OUTPUT_FOLDER <- if (length(args) >= 2) args[2] else "outputs/measurement"
BOOTSTRAPS <- 10000L
SEED <- 20260826L
 
# ------------------------ 1. PACKAGES AND FOLDERS ---------------------------
needed <- c("lavaan", "semTools", "semPlot", "readr", "psych")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Install these packages once, then rerun: install.packages(c(",
    paste(sprintf("'%s'", missing), collapse = ", "), "))"
  )
}
 
suppressPackageStartupMessages({
  library(lavaan)
  library(semTools)
  library(semPlot)
  library(readr)
  library(psych)
})
 
if (!file.exists(DATA_PATH)) stop(paste("Dataset not found:", DATA_PATH))
 
OUT <- normalizePath(OUTPUT_FOLDER, winslash = "/", mustWork = FALSE)
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
 
d_raw <- as.data.frame(readr::read_csv(DATA_PATH, show_col_types = FALSE),
                       check.names = FALSE)
 
# ------------------------------ 2. HELPERS ---------------------------------
require_cols <- function(data, variables, label = "variables") {
  absent <- setdiff(variables, names(data))
  if (length(absent) > 0) {
    stop(paste("Missing", label, ":", paste(absent, collapse = ", ")))
  }
}
 
fit_row <- function(fit, model_name, status = "Estimated") {
  fm <- fitMeasures(fit)
  estimator <- toupper(lavInspect(fit, "options")$estimator)
  if (estimator %in% c("WLSMV", "DWLS", "WLSM", "WLSMVS")) {
    variant <- "Scaled WLSMV/DWLS"
    candidates <- list(Chi_square=c("chisq.scaled","chisq"),df=c("df.scaled","df"),p=c("pvalue.scaled","pvalue"),CFI=c("cfi.scaled","cfi"),TLI=c("tli.scaled","tli"),RMSEA=c("rmsea.scaled","rmsea"),RMSEA_CI_low=c("rmsea.ci.lower.scaled","rmsea.ci.lower"),RMSEA_CI_high=c("rmsea.ci.upper.scaled","rmsea.ci.upper"),SRMR=c("srmr"))
  } else if (estimator %in% c("MLR", "MLM", "MLMV")) {
    variant <- "Robust ML"
    candidates <- list(Chi_square=c("chisq.scaled","chisq"),df=c("df.scaled","df"),p=c("pvalue.scaled","pvalue"),CFI=c("cfi.robust","cfi.scaled","cfi"),TLI=c("tli.robust","tli.scaled","tli"),RMSEA=c("rmsea.robust","rmsea.scaled","rmsea"),RMSEA_CI_low=c("rmsea.ci.lower.robust","rmsea.ci.lower.scaled","rmsea.ci.lower"),RMSEA_CI_high=c("rmsea.ci.upper.robust","rmsea.ci.upper.scaled","rmsea.ci.upper"),SRMR=c("srmr"))
  } else {
    variant <- "Standard ML"
    candidates <- list(Chi_square=c("chisq"),df=c("df"),p=c("pvalue"),CFI=c("cfi"),TLI=c("tli"),RMSEA=c("rmsea"),RMSEA_CI_low=c("rmsea.ci.lower"),RMSEA_CI_high=c("rmsea.ci.upper"),SRMR=c("srmr"))
  }
  pick <- function(x) { hit <- x[x %in% names(fm)]; if(length(hit)) unname(fm[hit[1]]) else NA_real_ }
  data.frame(Model=model_name,Status=status,Estimator=estimator,Fit_index_variant=variant,Chi_square=pick(candidates$Chi_square),df=pick(candidates$df),p=pick(candidates$p),CFI=pick(candidates$CFI),TLI=pick(candidates$TLI),RMSEA=pick(candidates$RMSEA),RMSEA_CI_low=pick(candidates$RMSEA_CI_low),RMSEA_CI_high=pick(candidates$RMSEA_CI_high),SRMR=pick(candidates$SRMR),row.names=NULL)
}
 
fit_or_status <- function(model_syntax, data, model_name,
                          estimator = "MLR", bootstrap = FALSE) {
  result <- tryCatch({
    if (bootstrap) {
      sem(model_syntax, data = data, estimator = "ML",
          se = "bootstrap", bootstrap = BOOTSTRAPS, fixed.x = TRUE,
          missing = "fiml")
    } else {
      sem(model_syntax, data = data, estimator = estimator, fixed.x = TRUE,
          missing = "fiml")
    }
  }, error = function(e) e)
 
  if (inherits(result, "error")) {
    return(list(
      fit = NULL,
      row = data.frame(
        Model = model_name,
        Status = paste("Not estimated:", result$message),
        Estimator = estimator, Fit_index_variant = NA_character_,
        Chi_square = NA, df = NA, p = NA, CFI = NA, TLI = NA,
        RMSEA = NA, RMSEA_CI_low = NA, RMSEA_CI_high = NA, SRMR = NA
      )
    ))
  }
  list(fit = result, row = fit_row(result, model_name))
}
 
save_sem_estimates <- function(fit, filename) {
  if (is.null(fit)) return(invisible(NULL))
  pe <- parameterEstimates(
    fit, standardized = TRUE, ci = TRUE,
    boot.ci.type = "bca.simple"
  )
  keep <- pe[pe$op %in% c("=~", "~", ":="),
             c("lhs", "op", "rhs", "est", "se", "ci.lower", "ci.upper",
               "pvalue", "std.all")]
  names(keep) <- c(
    "Left_hand_side", "Operator", "Right_hand_side", "Estimate", "SE",
    "CI_low", "CI_high", "p", "Standardized_estimate"
  )
  write.csv(keep, filename, row.names = FALSE)
}
 
write_png <- function(path, plot_fun, width = 14, height = 10) {
  png(path, width = width, height = height, units = "in", res = 600)
  tryCatch(plot_fun(), finally = dev.off())
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
    stop(paste("Figure was not written correctly:", path))
  }
}
 
plot_loading_profile <- function(tab, construct_col, item_col, loading_col,
                                 title, path) {
  z <- tab[, c(construct_col, item_col, loading_col)]
  names(z) <- c("Construct", "Item", "Loading")
  z <- z[is.finite(z$Loading), , drop = FALSE]
  z <- z[order(z$Construct, z$Loading), , drop = FALSE]
  cols <- setNames(grDevices::hcl.colors(length(unique(z$Construct)), "Dark 3"),
                   unique(z$Construct))
  write_png(path, function() {
    old <- par(mar = c(5, 11, 4, 2)); on.exit(par(old), add = TRUE)
    y <- seq_len(nrow(z))
    plot(z$Loading, y, xlim = c(0, 1), yaxt = "n", ylab = "",
         xlab = "Standardized factor loading", pch = 19,
         col = cols[z$Construct], main = title)
    axis(2, at = y, labels = z$Item, las = 1, cex.axis = .72)
    abline(v = c(.50, .70), lty = c(3, 2), col = c("grey55", "grey30"))
    legend("bottomright", legend = names(cols), col = cols, pch = 19,
           bty = "n", cex = .78)
    grid(nx = NA, ny = NULL, col = "grey92")
  }, width = 11, height = max(8, .22 * nrow(z) + 2.5))
}
 
plot_fit_comparison <- function(tab, title, path) {
  metrics <- c("CFI", "TLI", "RMSEA", "SRMR")
  vals <- as.matrix(tab[, metrics, drop = FALSE])
  rownames(vals) <- tab$Model
  write_png(path, function() {
    old <- par(mfrow = c(2, 2), mar = c(7, 4, 3, 1)); on.exit(par(old), add = TRUE)
    for (j in seq_along(metrics)) {
      lim <- if (metrics[j] %in% c("CFI", "TLI")) c(0, 1) else c(0, max(.12, vals[, j], na.rm = TRUE) * 1.10)
      barplot(vals[, j], names.arg = rownames(vals), las = 2, ylim = lim,
              col = "#52796F", border = NA, ylab = metrics[j], main = metrics[j])
      if (metrics[j] %in% c("CFI", "TLI")) abline(h = .90, lty = 3, col = "grey35")
      if (metrics[j] == "RMSEA") abline(h = .08, lty = 3, col = "grey35")
      if (metrics[j] == "SRMR") abline(h = .08, lty = 3, col = "grey35")
    }
    mtext(title, outer = TRUE, line = -1.2, cex = 1.15, font = 2)
  }, width = 12, height = 9)
}
 
# Classical alpha is reported alongside CFA-based reliability. Alpha is not
# substituted for CR or AVE; the three quantities answer different questions.
alpha_safe <- function(data, items) {
  require_cols(data, items, "alpha items")
  x <- data[, items, drop = FALSE]
  x[] <- lapply(x, as.numeric)
  suppressWarnings(as.numeric(psych::alpha(x, check.keys = FALSE,
                                             warnings = FALSE)$total$raw_alpha))
}
 
# Reliability from completely standardized loadings. For standardized
# solutions, 1-lambda^2 is the indicator residual variance used here.
reliability_from_fit <- function(fit, model_name, alpha_lookup = list()) {
  if (is.null(fit)) return(data.frame())
  ss <- standardizedSolution(fit)
  sl <- ss[ss$op == "=~" & is.finite(ss$est.std),
           c("lhs", "rhs", "est.std"), drop = FALSE]
  if (!nrow(sl)) return(data.frame())
  do.call(rbind, lapply(unique(sl$lhs), function(f) {
    z <- sl$est.std[sl$lhs == f]
    av <- if (f %in% names(alpha_lookup)) alpha_lookup[[f]] else NA_real_
    data.frame(
      Model = model_name, Construct = f, Indicators = length(z),
      Cronbach_alpha = av,
      Composite_reliability = sum(z)^2 / (sum(z)^2 + sum(1 - z^2)),
      AVE = mean(z^2), Loading_min = min(z), Loading_max = max(z),
      Alpha_note = ifelse(is.finite(av),
                          "Raw-item alpha",
                          "Not applicable/not calculated for latent higher-order or parcel specification"),
      row.names = NULL
    )
  }))
}
 
extract_loadings <- function(fit, model_name) {
  if (is.null(fit)) return(data.frame())
  pe <- parameterEstimates(fit, standardized = TRUE, ci = TRUE)
  pe <- pe[pe$op == "=~",
           c("lhs", "rhs", "est", "se", "z", "pvalue", "ci.lower",
             "ci.upper", "std.all"), drop = FALSE]
  names(pe) <- c("Construct", "Indicator", "Estimate", "SE", "z", "p",
                 "CI_low", "CI_high", "Standardized_loading")
  pe$Model <- model_name
  pe[, c("Model", setdiff(names(pe), "Model")), drop = FALSE]
}
 
model_admissibility <- function(fit, model_name) {
  if (is.null(fit)) {
    return(data.frame(Model=model_name, Converged=FALSE,
      Post_check_passed=FALSE, Negative_residual_variances=NA,
      Standardized_loadings_abs_gt_1=NA, R2_outside_0_1=NA))
  }
  pe <- parameterEstimates(fit, standardized = TRUE)
  neg_resid <- sum(pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0,
                   na.rm = TRUE)
  bad_loading <- sum(pe$op == "=~" & abs(pe$std.all) > 1, na.rm = TRUE)
  r2 <- tryCatch(unname(lavInspect(fit, "r2")), error = function(e) numeric(0))
  data.frame(
    Model = model_name,
    Converged = isTRUE(lavInspect(fit, "converged")),
    Post_check_passed = isTRUE(lavInspect(fit, "post.check")),
    Negative_residual_variances = neg_resid,
    Standardized_loadings_abs_gt_1 = bad_loading,
    R2_outside_0_1 = sum(r2 < 0 | r2 > 1, na.rm = TRUE),
    row.names = NULL
  )
}
 
model_r2 <- function(fit, model_name) {
  if (is.null(fit)) return(data.frame())
  z <- tryCatch(lavInspect(fit, "r2"), error = function(e) numeric(0))
  if (!length(z)) return(data.frame())
  data.frame(Model=model_name, Endogenous_variable=names(z), R_squared=unname(z),
             row.names=NULL)
}
 
# Create three balanced parcels. For a five-item scale, the parcels contain
# items (1,4), (2,5), and (3), preventing a single parcel from carrying all
# early or late items in the scale.
make_three_parcels <- function(data, item_names, prefix) {
  require_cols(data, item_names, prefix)
  groups <- list(item_names[c(1, 4)], item_names[c(2, 5)], item_names[3])
  if (length(item_names) == 7) {
    groups <- list(item_names[c(1, 4, 7)], item_names[c(2, 5)], item_names[c(3, 6)])
  }
  if (length(item_names) == 8) {
    groups <- list(item_names[c(1, 4, 7)], item_names[c(2, 5, 8)], item_names[c(3, 6)])
  }
  for (i in 1:3) {
    data[[paste0(prefix, "_p", i)]] <- rowMeans(data[, groups[[i]], drop = FALSE], na.rm = FALSE)
  }
  data
}
 
# ------------------------------ 3. VARIABLES --------------------------------
ccbs <- paste0("CCB_", 1:5)
lmx1 <- paste0("LMX1_", 1:7)
lmx2 <- paste0("LMX2_", 1:7)
task <- paste0("TP_", 1:5)
contextual <- paste0("CP_", 1:8)
cwb <- paste0("CWB_", 1:5)
ordinal_items <- c(ccbs, lmx1, lmx2, task, contextual, cwb)
require_cols(d_raw, ordinal_items, "measurement items")
 
# ==========================================================================
# PART A. ORDINAL MEASUREMENT MODEL (WLSMV)
# ==========================================================================
d_ordinal <- d_raw
for (v in ordinal_items) d_ordinal[[v]] <- ordered(d_ordinal[[v]])
 
six_factor_cfa <- paste(c(
  paste("CCB =~", paste(ccbs, collapse = " + ")),
  paste("LMX_T1 =~", paste(lmx1, collapse = " + ")),
  paste("LMX_T2 =~", paste(lmx2, collapse = " + ")),
  paste("Task =~", paste(task, collapse = " + ")),
  paste("Contextual =~", paste(contextual, collapse = " + ")),
  paste("CWB =~", paste(cwb, collapse = " + "))
), collapse = "\n")
 
one_factor_cfa <- paste("General =~", paste(ordinal_items, collapse = " + "))
 
fit_six <- cfa(
  six_factor_cfa, data = d_ordinal, ordered = ordinal_items,
  estimator = "WLSMV", parameterization = "theta"
)
fit_one <- cfa(
  one_factor_cfa, data = d_ordinal, ordered = ordinal_items,
  estimator = "WLSMV", parameterization = "theta"
)
 
write.csv(
  rbind(
    fit_row(fit_six, "Six-factor ordinal CFA"),
    fit_row(fit_one, "One-factor ordinal CFA")
  ),
  file.path(TAB, "Table_R1_ordinal_CFA_fit.csv"), row.names = FALSE
)
 
loadings <- parameterEstimates(fit_six, standardized = TRUE)
loadings <- loadings[loadings$op == "=~", c("lhs", "rhs", "std.all", "pvalue")]
names(loadings) <- c("Construct", "Item", "Standardized_loading", "p")
write.csv(loadings, file.path(TAB, "Table_R2_ordinal_CFA_standardized_loadings.csv"), row.names = FALSE)
 
# Measurement figures. The loading profile is easier to read in the thesis
# than a single crowded diagram containing all 37 observed variables.
plot_loading_profile(
  loadings, "Construct", "Item", "Standardized_loading",
  "Six-factor ordinal CFA: standardized item loadings",
  file.path(FIG, "Figure_R5_ordinal_CFA_standardized_loadings.png")
)
cfa_fit_for_plot <- rbind(
  fit_row(fit_six, "Six-factor"),
  fit_row(fit_one, "One-factor")
)
plot_fit_comparison(
  cfa_fit_for_plot,
  "Ordinal CFA fit comparison",
  file.path(FIG, "Figure_R6_six_factor_vs_one_factor_CFA_fit.png")
)
 
# --------------------------------------------------------------------------
# LMX response-category audit (invariance is estimated once in Part C4)
# --------------------------------------------------------------------------
category_audit <- do.call(rbind, lapply(c(lmx1, lmx2), function(v) {
  tab <- table(d_ordinal[[v]], useNA = "no")
  data.frame(Item=v, Observed_categories=paste(names(tab),collapse=", "),
    Number_of_categories=length(tab), Number_of_thresholds=max(length(tab)-1,0),
    Minimum_category_n=ifelse(length(tab)>0,min(tab),NA), row.names=NULL)
}))
write.csv(category_audit,file.path(TAB,"Table_R3_LMX_category_audit.csv"),row.names=FALSE)
 
# ==========================================================================
# PART B. LATENT SEM WITH CONTINUOUS PARCELS (MLR + BOOTSTRAP)
# ==========================================================================
# The SEM section is intentionally separate from the ordinal CFA section.
# It uses balanced item parcels and conscientiousness facet scores as continuous
# indicators, so that the full structural model is estimable and can be drawn
# as an AMOS-style path diagram.
 
d_sem <- d_raw
numeric_items <- c(ordinal_items, grep("^(SE|OR|DU|AS|SD|CA)_", names(d_sem), value = TRUE))
for (v in numeric_items) d_sem[[v]] <- as.numeric(d_sem[[v]])
 
# Reverse-code each CWB item before making CWB-R parcels. Higher CWB-R means
# less counterproductive behaviour and therefore better work performance.
if (any(unlist(d_sem[cwb]) < 1 | unlist(d_sem[cwb]) > 5, na.rm = TRUE)) {
  stop("CWB items must use the 1–5 response range before reverse coding.")
}
for (i in seq_along(cwb)) d_sem[[paste0("CWBR_", i)]] <- 6 - d_sem[[cwb[i]]]
cwbr <- paste0("CWBR_", 1:5)
 
d_sem <- make_three_parcels(d_sem, ccbs, "CCB")
d_sem <- make_three_parcels(d_sem, lmx1, "LMX1")
d_sem <- make_three_parcels(d_sem, lmx2, "LMX2")
d_sem <- make_three_parcels(d_sem, task, "TASK")
d_sem <- make_three_parcels(d_sem, contextual, "CONTEXT")
d_sem <- make_three_parcels(d_sem, cwbr, "CWBR")
 
# Six conscientiousness facets are used as indicators of the latent trait.
facet_prefixes <- c("SE", "OR", "DU", "AS", "SD", "CA")
facet_scores <- character(0)
for (facet in facet_prefixes) {
  its <- grep(paste0("^", facet, "_"), names(d_sem), value = TRUE)
  if (length(its) < 2) stop(paste("At least two items are required for facet", facet))
  score_name <- paste0(facet, "_facet")
  d_sem[[score_name]] <- rowMeans(d_sem[, its, drop = FALSE], na.rm = FALSE)
  facet_scores <- c(facet_scores, score_name)
}
 
# Covariates. Education is coded through two indicators, treating the lowest
# observed education category as the reference category.
covariates <- character(0)
for (v in c("Gender_Female", "Marital_Status_Married")) {
  if (v %in% names(d_sem)) {
    d_sem[[v]] <- as.numeric(d_sem[[v]])
    covariates <- c(covariates, v)
  }
}
 
if ("Education_Level" %in% names(d_sem)) {
  ed <- sort(unique(d_sem$Education_Level[!is.na(d_sem$Education_Level)]))
  if (length(ed) >= 2) {
    for (j in 2:length(ed)) {
      nm <- paste0("Education_", j)
      d_sem[[nm]] <- as.numeric(d_sem$Education_Level == ed[j])
      d_sem[[nm]][is.na(d_sem$Education_Level)] <- NA
      covariates <- c(covariates, nm)
    }
  }
}
 
for (v in c("Total_Experience_Years", "Tenure_Org_Months", "Tenure_Supervisor_Months")) {
  if (v %in% names(d_sem)) {
    nm <- paste0("z_", v)
    d_sem[[nm]] <- as.numeric(scale(d_sem[[v]]))
    covariates <- c(covariates, nm)
  }
}
 
control_rhs <- if (length(covariates)) paste(covariates, collapse = " + ") else ""
with_controls <- function(base) if (nchar(control_rhs)) paste0(base, " + ", control_rhs) else base
 
# ----------------------------------------------------------------------------
# B1. First-order performance SEM
# ----------------------------------------------------------------------------
first_order_sem_model <- paste(c(
  "CCB =~ CCB_p1 + CCB_p2 + CCB_p3",
  "LMX_T1 =~ LMX1_p1 + LMX1_p2 + LMX1_p3",
  "LMX_T2 =~ LMX2_p1 + LMX2_p2 + LMX2_p3",
  "Conscientiousness =~ SE_facet + OR_facet + DU_facet + AS_facet + SD_facet + CA_facet",
  "Task =~ TASK_p1 + TASK_p2 + TASK_p3",
  "Contextual =~ CONTEXT_p1 + CONTEXT_p2 + CONTEXT_p3",
  "Reverse_CWB =~ CWBR_p1 + CWBR_p2 + CWBR_p3",
  with_controls("LMX_T2 ~ a*CCB + LMX_T1 + Conscientiousness"),
  with_controls("Task ~ b_task*LMX_T2 + c_task*CCB + Conscientiousness"),
  with_controls("Contextual ~ b_context*LMX_T2 + c_context*CCB + Conscientiousness"),
  with_controls("Reverse_CWB ~ b_cwb*LMX_T2 + c_cwb*CCB + Conscientiousness"),
  "indirect_task := a*b_task",
  "indirect_contextual := a*b_context",
  "indirect_reverse_CWB := a*b_cwb"
), collapse = "\n")
 
first_mlr <- fit_or_status(first_order_sem_model, d_sem,
                           "First-order performance-dimensions SEM: MLR")
first_boot <- fit_or_status(first_order_sem_model, d_sem,
                            "First-order performance-dimensions SEM: ML bootstrap",
                            bootstrap = TRUE)
write.csv(rbind(first_mlr$row, first_boot$row),
          file.path(TAB, "Table_R5_first_order_SEM_fit.csv"), row.names = FALSE)
save_sem_estimates(first_boot$fit, file.path(TAB, "Table_R6_first_order_SEM_estimates_bootstrap.csv"))
writeLines(first_order_sem_model, file.path(OUT, "first_order_SEM_model_syntax.txt"))
 
# ----------------------------------------------------------------------------
# B2. Second-order IWP SEM
# ----------------------------------------------------------------------------
second_order_sem_model <- paste(c(
  "CCB =~ CCB_p1 + CCB_p2 + CCB_p3",
  "LMX_T1 =~ LMX1_p1 + LMX1_p2 + LMX1_p3",
  "LMX_T2 =~ LMX2_p1 + LMX2_p2 + LMX2_p3",
  "Conscientiousness =~ SE_facet + OR_facet + DU_facet + AS_facet + SD_facet + CA_facet",
  "Task =~ TASK_p1 + TASK_p2 + TASK_p3",
  "Contextual =~ CONTEXT_p1 + CONTEXT_p2 + CONTEXT_p3",
  "Reverse_CWB =~ CWBR_p1 + CWBR_p2 + CWBR_p3",
  "IWP =~ Task + Contextual + Reverse_CWB",
  with_controls("LMX_T2 ~ a*CCB + LMX_T1 + Conscientiousness"),
  with_controls("IWP ~ b*LMX_T2 + cprime*CCB + Conscientiousness"),
  "indirect := a*b",
  "total := cprime + indirect"
), collapse = "\n")
 
second_mlr <- fit_or_status(second_order_sem_model, d_sem,
                            "Second-order IWP SEM: MLR")
second_boot <- fit_or_status(second_order_sem_model, d_sem,
                             "Second-order IWP SEM: ML bootstrap",
                             bootstrap = TRUE)
write.csv(rbind(second_mlr$row, second_boot$row),
          file.path(TAB, "Table_R7_second_order_IWP_SEM_fit.csv"), row.names = FALSE)
save_sem_estimates(second_boot$fit, file.path(TAB, "Table_R8_second_order_IWP_SEM_estimates_bootstrap.csv"))
writeLines(second_order_sem_model, file.path(OUT, "second_order_IWP_SEM_model_syntax.txt"))
 
# ----------------------------------------------------------------------------
# B3. AMOS-style figures
# ----------------------------------------------------------------------------
if (!is.null(first_mlr$fit)) {
  write_png(
    file.path(FIG, "Figure_R1_first_order_SEM_full_measurement_and_structural.png"),
    function() { semPaths(
      first_mlr$fit, what = "std", whatLabels = "std", style = "lisrel",
      layout = "tree2", rotation = 2, residuals = FALSE, intercepts = FALSE,
      exoCov = FALSE, nCharNodes = 0, edge.label.cex = .55,
      sizeLat = 8, sizeMan = 3.2, mar = c(5, 5, 5, 5)
    ); title("First-order SEM: measurement and structural model", line=3) },
    width = 18, height = 12
  )
  write_png(
    file.path(FIG, "Figure_R2_first_order_SEM_structural_paths.png"),
    function() { semPaths(
      first_mlr$fit, what = "std", whatLabels = "std", style = "lisrel",
      structural = TRUE, layout = "tree2", rotation = 2,
      residuals = FALSE, intercepts = FALSE, exoCov = FALSE,
      nCharNodes = 0, edge.label.cex = .90, sizeLat = 11, sizeMan = 6,
      mar = c(6, 6, 6, 6)
    ); title("First-order SEM: structural paths", line=3) },
    width = 14, height = 9
  )
}
 
if (!is.null(second_mlr$fit)) {
  write_png(
    file.path(FIG, "Figure_R3_second_order_IWP_SEM_full_measurement_and_structural.png"),
    function() { semPaths(
      second_mlr$fit, what = "std", whatLabels = "std", style = "lisrel",
      layout = "tree2", rotation = 2, residuals = FALSE, intercepts = FALSE,
      exoCov = FALSE, nCharNodes = 0, edge.label.cex = .55,
      sizeLat = 8, sizeMan = 3.2, mar = c(5, 5, 5, 5)
    ); title("Second-order IWP SEM: measurement and structural model", line=3) },
    width = 18, height = 12
  )
  write_png(
    file.path(FIG, "Figure_R4_second_order_IWP_SEM_structural_paths.png"),
    function() { semPaths(
      second_mlr$fit, what = "std", whatLabels = "std", style = "lisrel",
      structural = TRUE, layout = "tree2", rotation = 2,
      residuals = FALSE, intercepts = FALSE, exoCov = FALSE,
      nCharNodes = 0, edge.label.cex = .90, sizeLat = 11, sizeMan = 6,
      mar = c(6, 6, 6, 6)
    ); title("Second-order IWP SEM: structural paths", line=3) },
    width = 14, height = 9
  )
}
 
 
 
# ==========================================================================
# PART C. COMPLETION ANALYSES REQUIRED BY CHAPTER 3
# ==========================================================================
 
# C1. CFA-based CR, AVE, Fornell-Larcker and HTMT.
ss <- standardizedSolution(fit_six)
sl <- ss[ss$op == "=~", c("lhs", "rhs", "est.std")]
rel <- do.call(rbind, lapply(unique(sl$lhs), function(f) {
  z <- sl$est.std[sl$lhs == f]
  data.frame(Construct=f, Composite_reliability=sum(z)^2/(sum(z)^2+sum(1-z^2)),
             AVE=mean(z^2), Loading_min=min(z), Loading_max=max(z))
}))
write.csv(rel, file.path(TAB,"Table_R9_CFA_CR_AVE.csv"), row.names=FALSE)
fc <- lavInspect(fit_six,"cor.lv"); fl <- fc
diag(fl) <- sqrt(rel$AVE[match(rownames(fl),rel$Construct)])
write.csv(cbind(Construct=rownames(fl),fl),file.path(TAB,"Table_R10_Fornell_Larcker.csv"),row.names=FALSE)
Rpoly <- lavCor(d_ordinal[,ordinal_items],ordered=ordinal_items,output="cor")
blocks <- list(CCB=ccbs,LMX_T1=lmx1,LMX_T2=lmx2,Task=task,Contextual=contextual,CWB=cwb)
htmt <- function(a,b){h<-mean(abs(Rpoly[a,b])); ma<-mean(abs(Rpoly[a,a][upper.tri(Rpoly[a,a])])); mb<-mean(abs(Rpoly[b,b][upper.tri(Rpoly[b,b])])); h/sqrt(ma*mb)}
hm <- matrix(1,6,6,dimnames=list(names(blocks),names(blocks)))
for(i in 1:6) for(j in 1:6) if(j>i) hm[i,j]<-hm[j,i]<-htmt(blocks[[i]],blocks[[j]])
write.csv(cbind(Construct=rownames(hm),hm),file.path(TAB,"Table_R11_HTMT.csv"),row.names=FALSE)
 
# C2. KMO, Bartlett, Harman one-factor variance, and MLR CFA sensitivity.
k <- psych::KMO(Rpoly); bt <- psych::cortest.bartlett(Rpoly,n=nrow(d_raw)); h <- psych::principal(d_raw[,ordinal_items],nfactors=1,rotate="none")
write.csv(data.frame(Diagnostic=c("KMO","Bartlett_chisq","Bartlett_df","Bartlett_p","Harman_first_factor"),
 Value=c(k$MSA,bt$chisq,bt$df,bt$p,h$Vaccounted[2,1])),file.path(TAB,"Table_R12_factorability_common_method.csv"),row.names=FALSE)
d_mlr<-d_raw; for(v in ordinal_items) d_mlr[[v]]<-as.numeric(d_mlr[[v]])
fit_mlr_cfa<-cfa(six_factor_cfa,data=d_mlr,estimator="MLR",missing="fiml")
write.csv(rbind(fit_row(fit_six,"Ordinal WLSMV"),fit_row(fit_mlr_cfa,"MLR sensitivity")),file.path(TAB,"Table_R13_CFA_estimator_sensitivity.csv"),row.names=FALSE)
 
# C3. Six-facet conscientiousness ordinal CFA.
con_items <- unlist(lapply(c("SE","OR","DU","AS","SD","CA"),function(x) paste0(x,"_",1:10)))
require_cols(d_raw,con_items,"conscientiousness items"); dcon<-d_raw; for(v in con_items) dcon[[v]]<-ordered(dcon[[v]])
con_model<-paste(sapply(c("SE","OR","DU","AS","SD","CA"),function(x) paste(x,"=~",paste(paste0(x,"_",1:10),collapse=" + "))),collapse="\n")
fit_con<-cfa(con_model,data=dcon,ordered=con_items,estimator="WLSMV",parameterization="theta")
write.csv(fit_row(fit_con,"Six-facet conscientiousness CFA"),file.path(TAB,"Table_R14_conscientiousness_CFA_fit.csv"),row.names=FALSE)
cl<-parameterEstimates(fit_con,standardized=TRUE); cl<-cl[cl$op=="=~",c("lhs","rhs","std.all","pvalue")]
write.csv(cl,file.path(TAB,"Table_R15_conscientiousness_loadings.csv"),row.names=FALSE)
cl_plot <- cl
names(cl_plot) <- c("Construct", "Item", "Standardized_loading", "p")
plot_loading_profile(
  cl_plot, "Construct", "Item", "Standardized_loading",
  "Conscientiousness six-facet CFA: standardized item loadings",
  file.path(FIG, "Figure_R7_conscientiousness_CFA_loadings.png")
)
 
# C4. Final ordinal LMX invariance using semTools identification.
lmx_config <- paste(paste("LMX_T1 =~",paste(lmx1,collapse=" + ")),paste("LMX_T2 =~",paste(lmx2,collapse=" + ")),paste(paste(lmx1,"~~",lmx2),collapse="\n"),sep="\n")
lfn<-list(LMX=c("LMX_T1","LMX_T2"))
longfit<-function(eq){syn<-semTools::measEq.syntax(configural.model=lmx_config,data=d_ordinal,ordered=c(lmx1,lmx2),parameterization="theta",ID.fac="std.lv",ID.cat="Wu.Estabrook.2016",longFacNames=lfn,long.equal=eq); cfa(as.character(syn),data=d_ordinal,ordered=c(lmx1,lmx2),estimator="WLSMV",parameterization="theta")}
lf0<-longfit(character(0)); lf1<-longfit("loadings"); lf2<-longfit(c("loadings","thresholds"))
linv<-rbind(fit_row(lf0,"Configural"),fit_row(lf1,"Metric"),fit_row(lf2,"Scalar/threshold")); linv$Delta_CFI<-c(NA,diff(linv$CFI)); linv$Delta_RMSEA<-c(NA,diff(linv$RMSEA)); linv$Delta_SRMR<-c(NA,diff(linv$SRMR))
write.csv(linv,file.path(TAB,"Table_R4_LMX_ordinal_invariance_FINAL.csv"),row.names=FALSE)
 
# Save the standardized loadings from the final invariance model as the
# authoritative equal-loading solution and visualize model-fit stability.
linv_load <- parameterEstimates(lf2, standardized=TRUE)
linv_load <- linv_load[linv_load$op=="=~",c("lhs","rhs","est","se","pvalue","std.all")]
names(linv_load) <- c("Construct","Item","Estimate","SE","p","Standardized_loading")
write.csv(linv_load,file.path(TAB,"Table_R16_LMX_scalar_invariance_loadings.csv"),row.names=FALSE)
plot_fit_comparison(
  linv,
  "Longitudinal LMX measurement-invariance fit",
  file.path(FIG, "Figure_R8_LMX_measurement_invariance_fit.png")
)
plot_loading_profile(
  linv_load, "Construct", "Item", "Standardized_loading",
  "LMX scalar/threshold invariance model: standardized loadings",
  file.path(FIG, "Figure_R9_LMX_invariance_standardized_loadings.png")
)
 
# C5. SEM estimator and control-block sensitivity.
noctrl <- gsub(paste0(" \\+ ",control_rhs),"",second_order_sem_model)
fit_noctrl <- fit_or_status(noctrl,d_sem,"Second-order IWP SEM without controls",estimator="MLR")
write.csv(rbind(second_mlr$row,second_boot$row,fit_noctrl$row),file.path(TAB,"Table_R17_SEM_estimator_control_sensitivity.csv"),row.names=FALSE)
save_sem_estimates(fit_noctrl$fit,file.path(TAB,"Table_R18_SEM_without_controls_estimates.csv"))
 
# Bootstrap indirect-effect plot and SEM fit comparison.
if (!is.null(first_boot$fit) || !is.null(second_boot$fit)) {
  effect_rows <- list()
  if (!is.null(first_boot$fit)) {
    p1 <- parameterEstimates(first_boot$fit, ci=TRUE, boot.ci.type="bca.simple")
    p1 <- p1[p1$op==":=",c("lhs","est","ci.lower","ci.upper","pvalue")]
    p1$Model <- "First-order"
    effect_rows[[length(effect_rows)+1]] <- p1
  }
  if (!is.null(second_boot$fit)) {
    p2 <- parameterEstimates(second_boot$fit, ci=TRUE, boot.ci.type="bca.simple")
    p2 <- p2[p2$op==":=",c("lhs","est","ci.lower","ci.upper","pvalue")]
    p2$Model <- "Second-order"
    effect_rows[[length(effect_rows)+1]] <- p2
  }
  indirect_plot <- do.call(rbind,effect_rows)
  names(indirect_plot)[1:5] <- c("Effect","Estimate","CI_low","CI_high","p")
  write.csv(indirect_plot,file.path(TAB,"Table_R19_bootstrap_defined_effects.csv"),row.names=FALSE)
  write_png(file.path(FIG,"Figure_R10_bootstrap_indirect_effects.png"),function(){
    z <- indirect_plot[is.finite(indirect_plot$Estimate) & is.finite(indirect_plot$CI_low) & is.finite(indirect_plot$CI_high),]
    labs <- paste(z$Model,z$Effect,sep=": "); y <- rev(seq_len(nrow(z)))
    old <- par(mar=c(5,13,4,2)); on.exit(par(old),add=TRUE)
    xr <- range(c(z$CI_low,z$CI_high,0),finite=TRUE)
    plot(z$Estimate,y,xlim=xr,yaxt="n",ylab="",xlab="Bootstrap effect estimate",
         pch=19,col="#354F52",main=paste0("Bootstrap defined effects (",format(BOOTSTRAPS,big.mark=",")," resamples)"))
    segments(z$CI_low,y,z$CI_high,y,lwd=2,col="#354F52")
    axis(2,at=y,labels=labs,las=1,cex.axis=.78); abline(v=0,lty=2,col="grey35"); grid(nx=NA,ny=NULL,col="grey92")
  },width=12,height=max(7,.55*nrow(indirect_plot)+2))
}
 
sem_fit_plot <- rbind(first_mlr$row, first_boot$row, second_mlr$row,
                      second_boot$row, fit_noctrl$row)
sem_fit_plot <- sem_fit_plot[sem_fit_plot$Status=="Estimated",,drop=FALSE]
if (nrow(sem_fit_plot)>0) plot_fit_comparison(
  sem_fit_plot,
  "SEM estimator and control-block fit comparison",
  file.path(FIG,"Figure_R11_SEM_fit_comparison.png")
)
 
# C5b. Fit-plausibility audit. This flags unusually strong fit without changing any estimate.
all_fit_tables <- rbind(
  transform(fit_row(fit_six,"Six-factor ordinal CFA"), Source="CFA"),
  transform(fit_row(fit_one,"One-factor ordinal CFA"), Source="CFA"),
  transform(fit_row(fit_con,"Conscientiousness CFA"), Source="CFA"),
  transform(fit_row(lf0,"LMX configural"), Source="Invariance"),
  transform(fit_row(lf1,"LMX metric"), Source="Invariance"),
  transform(fit_row(lf2,"LMX scalar/threshold"), Source="Invariance")
)
all_fit_tables$Near_perfect_fit_flag <- with(all_fit_tables, !is.na(CFI) & !is.na(RMSEA) & CFI >= .999 & RMSEA <= .001)
all_fit_tables$Interpretation <- ifelse(all_fit_tables$Near_perfect_fit_flag,"Inspect identification, response categories, residual correlations, and estimator sensitivity; do not alter the values.","No near-perfect-fit flag under the prespecified rule.")
write.csv(all_fit_tables,file.path(TAB,"Table_R20_fit_plausibility_audit.csv"),row.names=FALSE)
 
# C5c. Consolidated reliability, validity, loading, R-squared and admissibility
# outputs. These tables prevent first-order and higher-order coefficients from
# being mixed together during thesis reporting.
alpha_lookup <- list(
  CCB=alpha_safe(d_raw,ccbs), LMX_T1=alpha_safe(d_raw,lmx1),
  LMX_T2=alpha_safe(d_raw,lmx2), Task=alpha_safe(d_raw,task),
  Contextual=alpha_safe(d_raw,contextual), CWB=alpha_safe(d_raw,cwb),
  Reverse_CWB=alpha_safe(d_raw,cwb),
  Conscientiousness=alpha_safe(d_sem,facet_scores),
  SE=alpha_safe(d_raw,paste0("SE_",1:10)),
  OR=alpha_safe(d_raw,paste0("OR_",1:10)),
  DU=alpha_safe(d_raw,paste0("DU_",1:10)),
  AS=alpha_safe(d_raw,paste0("AS_",1:10)),
  SD=alpha_safe(d_raw,paste0("SD_",1:10)),
  CA=alpha_safe(d_raw,paste0("CA_",1:10))
)
 
complete_reliability <- do.call(rbind, list(
  reliability_from_fit(fit_six,"Ordinal six-factor CFA",alpha_lookup),
  reliability_from_fit(fit_con,"Conscientiousness six-facet CFA",alpha_lookup),
  reliability_from_fit(first_mlr$fit,"First-order parcel SEM",alpha_lookup),
  reliability_from_fit(second_mlr$fit,"Second-order IWP parcel SEM",alpha_lookup)
))
complete_reliability$CR_acceptable_070 <- complete_reliability$Composite_reliability >= .70
complete_reliability$AVE_acceptable_050 <- complete_reliability$AVE >= .50
write.csv(complete_reliability,
          file.path(TAB,"Table_R21_complete_alpha_CR_AVE.csv"),row.names=FALSE)
 
first_order_loadings <- extract_loadings(first_mlr$fit,"First-order parcel SEM")
write.csv(first_order_loadings,
          file.path(TAB,"Table_R22_first_order_SEM_loadings.csv"),row.names=FALSE)
 
second_all_loadings <- extract_loadings(second_mlr$fit,"Second-order IWP parcel SEM")
second_order_loadings <- second_all_loadings[
  second_all_loadings$Construct == "IWP" &
    second_all_loadings$Indicator %in% c("Task","Contextual","Reverse_CWB"),,
  drop=FALSE]
write.csv(second_order_loadings,
          file.path(TAB,"Table_R23_second_order_IWP_loadings.csv"),row.names=FALSE)
 
admissibility <- do.call(rbind,list(
  model_admissibility(fit_six,"Ordinal six-factor CFA"),
  model_admissibility(fit_con,"Conscientiousness six-facet CFA"),
  model_admissibility(lf2,"LMX scalar/threshold invariance"),
  model_admissibility(first_mlr$fit,"First-order parcel SEM"),
  model_admissibility(first_boot$fit,"First-order bootstrap SEM"),
  model_admissibility(second_mlr$fit,"Second-order IWP parcel SEM"),
  model_admissibility(second_boot$fit,"Second-order IWP bootstrap SEM")
))
write.csv(admissibility,file.path(TAB,"Table_R24_model_admissibility.csv"),row.names=FALSE)
 
r_squared <- do.call(rbind,list(
  model_r2(first_mlr$fit,"First-order parcel SEM"),
  model_r2(second_mlr$fit,"Second-order IWP parcel SEM")
))
write.csv(r_squared,file.path(TAB,"Table_R25_SEM_R_squared.csv"),row.names=FALSE)
 
# Additional thesis figures. R12 is an AMOS-like CFA measurement diagram;
# R13 isolates the higher-order IWP loadings; R14 summarizes CR and AVE.
write_png(file.path(FIG,"Figure_R12_six_factor_CFA_AMOS_measurement_model.png"),function(){
  semPaths(fit_six,what="std",whatLabels="std",style="lisrel",layout="tree2",
           rotation=2,residuals=FALSE,intercepts=FALSE,exoCov=TRUE,
           nCharNodes=0,edge.label.cex=.48,sizeLat=8,sizeMan=2.7,
           mar=c(5,5,7,5))
  title("Six-factor ordinal CFA measurement model",line=3)
},width=20,height=14)
 
if(nrow(second_order_loadings)>0) {
  plot_loading_profile(second_order_loadings,"Construct","Indicator",
                       "Standardized_loading",
                       "Second-order IWP standardized loadings",
                       file.path(FIG,"Figure_R13_second_order_IWP_loadings.png"))
}
 
plot_rel <- complete_reliability[
  complete_reliability$Model %in% c("Ordinal six-factor CFA",
                                     "Conscientiousness six-facet CFA"),,
  drop=FALSE]
write_png(file.path(FIG,"Figure_R14_CR_AVE_summary.png"),function(){
  labs <- paste(plot_rel$Construct,plot_rel$Model,sep="\n")
  vals <- rbind(plot_rel$Composite_reliability,plot_rel$AVE)
  old <- par(mar=c(11,5,4,2)); on.exit(par(old),add=TRUE)
  barplot(vals,beside=TRUE,names.arg=labs,las=2,ylim=c(0,1),
          col=c("#355070","#B56576"),border=NA,ylab="Coefficient",
          main="CFA-based composite reliability and average variance extracted")
  abline(h=c(.70,.50),lty=c(2,3),col=c("#355070","#B56576"))
  legend("topright",legend=c("Composite reliability","AVE"),
         fill=c("#355070","#B56576"),bty="n")
},width=15,height=10)
 
# C6. Prevent silent figure failures or SEM-figure overwriting and record reproducibility.
fp<-file.path(FIG,c("Figure_R1_first_order_SEM_full_measurement_and_structural.png","Figure_R3_second_order_IWP_SEM_full_measurement_and_structural.png"))
if(all(file.exists(fp)) && unname(tools::md5sum(fp)[1])==unname(tools::md5sum(fp)[2])) stop("First-order and second-order SEM figures are identical; inspect semPaths output.")
expected_figures <- paste0("Figure_R",1:14,c(
  "_first_order_SEM_full_measurement_and_structural.png",
  "_first_order_SEM_structural_paths.png",
  "_second_order_IWP_SEM_full_measurement_and_structural.png",
  "_second_order_IWP_SEM_structural_paths.png",
  "_ordinal_CFA_standardized_loadings.png",
  "_six_factor_vs_one_factor_CFA_fit.png",
  "_conscientiousness_CFA_loadings.png",
  "_LMX_measurement_invariance_fit.png",
  "_LMX_invariance_standardized_loadings.png",
  "_bootstrap_indirect_effects.png",
  "_SEM_fit_comparison.png",
  "_six_factor_CFA_AMOS_measurement_model.png",
  "_second_order_IWP_loadings.png",
  "_CR_AVE_summary.png"
))
figure_paths <- file.path(FIG,expected_figures)
figure_audit <- data.frame(
  File=expected_figures,
  Exists=file.exists(figure_paths),
  Bytes=ifelse(file.exists(figure_paths),file.info(figure_paths)$size,NA),
  Non_empty=ifelse(file.exists(figure_paths),file.info(figure_paths)$size>0,FALSE)
)
write.csv(figure_audit,file.path(OUT,"figure_integrity_audit.csv"),row.names=FALSE)
if(any(!figure_audit$Non_empty)) warning("One or more expected figures were not created. See figure_integrity_audit.csv.")
 
table_files <- sort(list.files(TAB,pattern="\\.csv$",full.names=TRUE))
table_audit <- data.frame(File=basename(table_files),Bytes=file.info(table_files)$size,
                          Non_empty=file.info(table_files)$size>0)
write.csv(table_audit,file.path(OUT,"table_integrity_audit.csv"),row.names=FALSE)
writeLines(c(paste("Dataset:",DATA_PATH),paste("N:",nrow(d_raw)),paste("Seed:",SEED),paste("Bootstrap:",BOOTSTRAPS),"Primary IWP: Task + Contextual + reverse-scored CWB","Only CCB T1 -> LMX T2 controlling LMX T1 is baseline-adjusted temporal."),file.path(OUT,"analysis_manifest_complete.txt"))
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_sessionInfo_complete.txt"))
 
# A concise user-facing record of what each output contains.
readme <- c(
  "CCB two-phase R/lavaan analysis output",
  "",
  "Table_R1: ordinal six-factor versus one-factor CFA fit.",
  "Table_R2: standardised ordinal CFA item loadings.",
  "Table_R3: LMX response-category audit.",
  "Table_R4: ordinal longitudinal LMX measurement invariance.",
  "Table_R5/R6: first-order SEM fit and bootstrap estimates.",
  "Table_R7/R8: second-order IWP SEM fit and bootstrap estimates.",
  "Table_R9: CFA-based composite reliability and AVE.",
  "Table_R10/R11: Fornell-Larcker and HTMT discriminant validity.",
  "Table_R12/R13: factorability/common-method diagnostics and estimator sensitivity.",
  "Table_R14/R15: conscientiousness CFA fit and item loadings.",
  "Table_R16: final LMX scalar/threshold invariance standardized loadings.",
  "Table_R17/R18: SEM estimator/control sensitivity and no-control estimates.",
  "Table_R19: bootstrap indirect and other defined effects.",
  "Table_R20: transparent near-perfect-fit plausibility flags; no fit value is modified.",
  "Table_R21: consolidated raw-item alpha, CFA/SEM composite reliability, and AVE.",
  "Table_R22: explicit first-order SEM standardized loadings.",
  "Table_R23: explicit higher-order Task, Contextual, and Reverse-CWB loadings on IWP.",
  "Table_R24: convergence, post-check, Heywood/loading, and R-squared admissibility audit.",
  "Table_R25: explained variance for SEM endogenous variables.",
  "",
  "Figures R1/R2: first-order SEM full and structural-only diagrams.",
  "Figures R3/R4: second-order IWP SEM full and structural-only diagrams.",
  "Figure R5: six-factor ordinal CFA standardized loading profile.",
  "Figure R6: six-factor versus one-factor CFA fit comparison.",
  "Figure R7: conscientiousness six-facet CFA loading profile.",
  "Figures R8/R9: LMX invariance fit and standardized-loading profile.",
  "Figure R10: bootstrap indirect and defined effects with confidence intervals.",
  "Figure R11: SEM estimator and control-block fit comparison.",
  "Figure R12: AMOS-like six-factor CFA measurement structure.",
  "Figure R13: isolated second-order IWP standardized loadings.",
  "Figure R14: CFA composite reliability and AVE summary.",
  "",
  "Interpretation rule: only CCB T1 -> LMX T2 is baseline-adjusted temporal.",
  "LMX T2 -> performance/IWP is a Phase-2 associational path.",
  "CWB was reverse-coded at the item level for the second-order IWP model."
)
writeLines(readme, file.path(OUT, "README_output_guide.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "R_sessionInfo.txt"))
 
cat("Completed. Results are saved in:\n", OUT, "\n")
