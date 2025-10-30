##################################################################
##              Comparing RF+ and RF- on lymphiods              ##
##################################################################
set.seed(123)

rm(list=ls())

modules <- readRDS("/Users/javad/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/Data/significant_modules_members.RDS")
suppressPackageStartupMessages({
  library(caret)
  library(dplyr)
  library(tibble)
  library(pROC)
  library(PRROC)
  library(glmnet)
  library(tidyr)
  library(ggplot2)
})



setwd("/Users/javad/OneDrive - University of Pittsburgh/RA/E-MTAB-6141(PEAC)/rdata/")

load("E-MTAB-6141.rdata")

mat <- readRDS("mat.RDS")

meta <- read.delim("metadata.tsv",sep="\t")

E.MTAB.6141.sdrf <- read.delim("~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/E-MTAB-6141(PEAC)//rdata/E-MTAB-6141.sdrf.txt")

table(E.MTAB.6141.sdrf$Characteristics.organism.part.)

intersect_ids <- intersect(colnames(mat),E.MTAB.6141.sdrf[,"Comment.GenentechID."])

sdrf_unique <- E.MTAB.6141.sdrf[!duplicated(E.MTAB.6141.sdrf[,"Comment.GenentechID."]),]

sdrf_common <- sdrf_unique[sdrf_unique[,"Comment.GenentechID."] %in%intersect_ids, ]


## =======================
## USER CONFIG
## =======================
USE_CLASS_WEIGHTS <- FALSE          # TRUE to enable glmnet 'weights'
SAMPLING_METHOD   <- "none"         # "none", "down", "up", or "smote"
K_TARGET          <- 10             # target number of folds
PATHOTYPE_FILTER  <- "lymphoid"     # set to NULL to use all pathotypes

## =======================
## Prep RF labels & alignment
## =======================
id_col <- "Comment.GenentechID."
mods_idx <- c(1, 5, 2, 10)
has_rf <- !is.na(sdrf_common$Characteristics.RF.)
sdrf_rf <- sdrf_common[has_rf, , drop = FALSE]
sdrf_rf$rf_group <- ifelse(sdrf_rf$Characteristics.RF. == 1, "RF+", "RF-")

if (!is.null(PATHOTYPE_FILTER)) {
  sdrf_rf <- sdrf_rf[sdrf_rf$Factor.Value.pathotype. == PATHOTYPE_FILTER, , drop = FALSE]
}

idx  <- match(colnames(mat), sdrf_rf[[id_col]])
keep <- !is.na(idx)
mat_rf        <- mat[, keep, drop = FALSE]
sdrf_rf_align <- sdrf_rf[idx[keep], , drop = FALSE]
stopifnot(identical(colnames(mat_rf), sdrf_rf_align[[id_col]]))
message("The mat_rf and sdrf_rf_align are aligned.")

# Labels: POSITIVE = RF+ (put positive first for caret)
# Encode RF+ as 1L, RF- as 0L so that factor(1,0) -> c("case","control") has RF+ as 'case'
y_all_num <- ifelse(sdrf_rf_align$rf_group == "RF+", 1L, 0L)
X_all <- t(mat_rf)  # samples x genes

## =======================
## Helpers
## =======================
# Build stratified folds with both classes in each holdout
make_safe_folds <- function(y, k = 5, max_attempts = 1000) {
  y <- factor(y)
  k <- max(2L, min(k, length(y)))
  for (i in seq_len(max_attempts)) {
    te_list <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
    ok <- all(vapply(te_list, function(te) nlevels(droplevels(y[te])) == 2L, logical(1)))
    if (ok) {
      return(list(
        index    = lapply(te_list, function(te) setdiff(seq_along(y), te)),
        indexOut = te_list
      ))
    }
  }
  stop("Couldn't make folds with both classes in each holdout. Try smaller k.")
}

# Custom summary: AUPRC, ROC, Sens, Spec, Precision, Recall, F1, Threshold
summary_pr_roc_f1 <- function(data, lev = NULL, model = NULL) {
  stopifnot(length(lev) == 2)
  pos <- lev[1]
  if (!pos %in% colnames(data)) stop("Positive class probability column not found in 'data'.")
  scores <- data[[pos]]
  y <- as.integer(data$obs == pos)
  
  roc_obj <- pROC::roc(response = y, predictor = scores, quiet = TRUE, direction = "<")
  auc_roc <- as.numeric(pROC::auc(roc_obj))
  
  pr <- PRROC::pr.curve(
    scores.class0 = scores[y == 1],  # positives
    scores.class1 = scores[y == 0],  # negatives
    curve = TRUE
  )
  auc_pr <- as.numeric(pr$auc.integral)
  
  ths <- unique(c(0, 0.01, 0.05, 0.1, pr$curve[,1], 0.9, 0.95, 0.99, 1))
  ths <- ths[ths >= 0 & ths <= 1]
  best <- list(f1 = -Inf, th = 0.5, sens = NA, spec = NA, ppv = NA, rec = NA)
  
  for (t in ths) {
    pred_pos <- scores >= t
    TP <- sum(pred_pos & y == 1)
    FP <- sum(pred_pos & y == 0)
    FN <- sum(!pred_pos & y == 1)
    TN <- sum(!pred_pos & y == 0)
    
    sens <- if ((TP+FN) > 0) TP/(TP+FN) else NA_real_     # recall
    ppv  <- if ((TP+FP) > 0) TP/(TP+FP) else NA_real_     # precision
    spec <- if ((TN+FP) > 0) TN/(TN+FP) else NA_real_
    f1   <- if (!is.na(sens) && !is.na(ppv) && (sens+ppv) > 0) 2*sens*ppv/(sens+ppv) else NA_real_
    
    if (!is.na(f1) && f1 > best$f1) best <- list(f1=f1, th=t, sens=sens, spec=spec, ppv=ppv, rec=sens)
  }
  
  if (!"pred" %in% names(data)) {
    neg <- setdiff(lev, pos)
    pred <- ifelse(scores >= best$th, pos, neg)
    data$pred <- factor(pred, levels = lev)
  }
  
  c(
    AUPRC     = auc_pr,
    ROC       = auc_roc,
    Sens      = best$sens,
    Spec      = best$spec,
    Precision = best$ppv,
    Recall    = best$rec,
    F1        = best$f1,
    Thr       = best$th
  )
}

## =======================
## Tuning grid
## =======================
grid <- expand.grid(
  alpha  = seq(0, 1, by = 0.25),
  lambda = 10^seq(-4, 1, length.out = 30)
)

## =======================
## Train per module (CV, stratified & safe)
## =======================
report <- tibble()
fits   <- list()

for (mod in mods_idx) {
  genes <- intersect(colnames(X_all), modules[[mod]])
  if (!length(genes)) next
  
  X_module <- X_all[, genes, drop = FALSE]
  
  nzv <- nearZeroVar(X_module)
  if (length(nzv)) X_module <- X_module[, -nzv, drop = FALSE]
  if (ncol(X_module) == 0L) next
  
  keep_rows <- complete.cases(X_module)
  Xm  <- X_module[keep_rows, , drop = FALSE]
  ymn <- y_all_num[keep_rows]
  
  # Factor outcome: POSITIVE FIRST -> "case" = RF+
  ym <- factor(ymn, levels = c(1L, 0L), labels = c("case","control"))
  if (nlevels(droplevels(ym)) < 2L) {
    report <- bind_rows(report, tibble(
      module       = mod,
      n_genes      = length(genes),
      n_after_nzv  = ncol(X_module),
      samples      = length(ym),
      pos_n        = sum(ym == "case"),
      neg_n        = sum(ym == "control"),
      pos_pct      = if (length(ym)) mean(ym == "case") else NA_real_,
      alpha = NA, lambda = NA,
      ROC = NA_real_, AUPRC = NA_real_, Sens = NA_real_, Spec = NA_real_,
      Precision = NA_real_, Recall = NA_real_, F1 = NA_real_, Thr = NA_real_,
      AUPRC_vs_prev = NA_character_
    ))
    next
  }
  
  # pick k safely based on class counts
  tab <- table(ym)
  k_here <- min(K_TARGET, as.integer(tab["case"]), as.integer(tab["control"]))
  if (is.na(k_here) || k_here < 2L) k_here <- 2L
  sf <- make_safe_folds(ym, k = k_here)
  
  # Optional class weights
  wts <- NULL
  if (USE_CLASS_WEIGHTS) {
    p_case <- mean(ym == "case")
    p_ctrl <- mean(ym == "control")
    wts <- ifelse(ym == "case", 0.5 / p_case, 0.5 / p_ctrl)
  }
  
  # Optional sampling
  samp <- if (SAMPLING_METHOD %in% c("down","up","smote")) SAMPLING_METHOD else NULL
  
  ctrl <- trainControl(
    method = "cv",
    number = length(sf$index),
    classProbs = TRUE,
    summaryFunction = summary_pr_roc_f1,
    savePredictions = "final",
    index = sf$index,
    indexOut = sf$indexOut,
    sampling = samp
  )
  
  fit <- train(
    x = as.matrix(Xm), y = ym,
    method = "glmnet",
    trControl = ctrl,
    tuneGrid  = grid,
    preProcess = c("center","scale"),
    metric = "AUPRC",
    weights = wts
  )
  
  best <- fit$results %>%
    dplyr::filter(alpha == fit$bestTune$alpha, lambda == fit$bestTune$lambda)
  
  pos_n   <- sum(ym == "case")
  neg_n   <- sum(ym == "control")
  pos_pct <- pos_n / (pos_n + neg_n)
  
  report <- bind_rows(report, tibble(
    module       = mod,
    n_genes      = length(genes),
    n_after_nzv  = ncol(X_module),
    samples      = length(ym),
    pos_n        = pos_n,
    neg_n        = neg_n,
    pos_pct      = pos_pct,
    alpha        = fit$bestTune$alpha,
    lambda       = fit$bestTune$lambda,
    ROC          = best$ROC,
    AUPRC        = best$AUPRC,
    Sens         = best$Sens,
    Spec         = best$Spec,
    Precision    = best$Precision,
    Recall       = best$Recall,
    F1           = best$F1,
    Thr          = best$Thr,
    AUPRC_vs_prev = sprintf("%.3f (π=%.1f%%)", best$AUPRC, 100*pos_pct)
  ))
  
  fits[[as.character(mod)]] <- fit
}

## =======================
## Final report
## =======================
report <- report %>%
  arrange(desc(AUPRC), desc(ROC)) %>%
  mutate(pos_pct = round(100*pos_pct, 1))

print(report)

## =======================
## Plot ROC & AUPRC per module
## =======================
df_long <- report %>%
  dplyr::select(module, ROC, AUPRC, F1) %>%
  tidyr::pivot_longer(cols = c("ROC", "AUPRC"),
                      names_to = "Metric",
                      values_to = "Value")

metric_colors <- c("ROC" = "#1f77b4", "AUPRC" = "#ff7f0e")

df_long$Value <- as.numeric(df_long$Value)
pos_rate <- mean(report$pos_n / (report$pos_n + report$neg_n))

write.table(df_long,file="~/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/E-MTAB-6141(PEAC)/Analysis/Permutation_RF+_vs_RF-.csv",
            col.names = F,
            row.names = F)

ggplot(df_long, aes(x = factor(module), y = Value, fill = Metric)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.7),
           width = 0.6, color = "black", size = 1) +
  geom_hline(yintercept = pos_rate, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = Inf, y = pos_rate,
           label = sprintf("Baseline = %.2f", pos_rate),
           color = "red", fontface = "bold", hjust = 1.02, vjust = -0.5, size = 4) +
  scale_fill_manual(values = metric_colors) +
  coord_cartesian(ylim = c(0.5, 1)) +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(),
        axis.title  = element_text(face = "bold", size = 14),
        axis.text   = element_text(face = "bold", size = 12)) +
  labs(x = "Module", y = "Performance Value",
       title = sprintf("Module-wise ROC & AUPRC (RF+ vs RF-%s)",
                       ifelse(is.null(PATHOTYPE_FILTER), "", paste0(", ", PATHOTYPE_FILTER))))




#################################################################
##                   Making permutation plot                   ##
#################################################################
suppressPackageStartupMessages({
  library(purrr); library(tidyr); library(ggplot2)
})

## =======================
## Permutation test setup
## =======================
B <- 200   # number of label shuffles per module (raise to ~1000 if time allows)

# helper to (re)fit a glmnet model and return metrics
fit_metrics_glmnet <- function(Xm, ym, grid, use_wts=USE_CLASS_WEIGHTS, samp=SAMPLING_METHOD, k_target=K_TARGET) {
  # pick k safely for current labels
  tab <- table(ym)
  k_here <- min(k_target, as.integer(tab["case"]), as.integer(tab["control"]))
  if (is.na(k_here) || k_here < 2L) k_here <- 2L
  sf <- make_safe_folds(ym, k = k_here)
  
  wts <- NULL
  if (use_wts) {
    p_case <- mean(ym == "case"); p_ctrl <- 1 - p_case
    wts <- ifelse(ym == "case", 0.5/p_case, 0.5/p_ctrl)
  }
  samp <- if (samp %in% c("down","up","smote")) samp else NULL
  
  ctrl <- trainControl(
    method="cv", number=length(sf$index),
    classProbs=TRUE, summaryFunction=summary_pr_roc_f1,
    savePredictions="final", index=sf$index, indexOut=sf$indexOut,
    sampling=samp
  )
  
  fit <- train(
    x=as.matrix(Xm), y=ym,
    method="glmnet", trControl=ctrl, tuneGrid=grid,
    preProcess=c("center","scale"),
    metric="AUPRC", weights=wts
  )
  
  best <- dplyr::filter(fit$results,
                        alpha==fit$bestTune$alpha,
                        lambda==fit$bestTune$lambda)
  c(AUPRC = best$AUPRC, ROC = best$ROC, F1 = best$F1)
}

# build module matrix + labels once
get_module_xy <- function(mod) {
  genes <- intersect(colnames(X_all), modules[[mod]])
  if (!length(genes)) return(NULL)
  X_module <- X_all[, genes, drop=FALSE]
  nzv <- nearZeroVar(X_module)
  if (length(nzv)) X_module <- X_module[, -nzv, drop=FALSE]
  if (ncol(X_module) == 0L) return(NULL)
  
  keep_rows <- complete.cases(X_module)
  Xm <- X_module[keep_rows, , drop=FALSE]
  ym <- factor(y_all_num[keep_rows], levels=c(1L,0L), labels=c("case","control"))
  list(Xm=Xm, ym=ym)
}

## =========================================================
## Run permutations per module and gather distributions
## =========================================================
mods_to_do <- report$module   # only modules that produced an observed result

perm_list <- list()
pval_list <- list()

for (mod in mods_to_do) {
  xy <- get_module_xy(mod)
  if (is.null(xy)) next
  
  Xm <- xy$Xm; ym_obs <- xy$ym
  
  # observed (use already-computed numbers for speed)
  obs_row <- dplyr::filter(report, module == mod)
  obs <- c(AUPRC = obs_row$AUPRC, ROC = obs_row$ROC)
  
  # permutations
  set.seed(1000 + as.integer(mod))
  perm_vals <- replicate(B, {
    y_perm <- sample(ym_obs)                # shuffle labels
    out <- try(fit_metrics_glmnet(Xm, y_perm, grid), silent=TRUE)
    if (inherits(out, "try-error")) c(AUPRC=NA_real_, ROC=NA_real_, F1=NA_real_) else out
  })
  
  perm_df <- as.data.frame(t(perm_vals))
  perm_df$module <- mod
  perm_list[[as.character(mod)]] <- perm_df
  
  # p-values (one-sided: larger is better)
  p_aucpr <- (1 + sum(perm_df$AUPRC >= obs["AUPRC"], na.rm=TRUE)) / (1 + sum(!is.na(perm_df$AUPRC)))
  p_roc   <- (1 + sum(perm_df$ROC   >= obs["ROC"],   na.rm=TRUE)) / (1 + sum(!is.na(perm_df$ROC)))
  pval_list[[as.character(mod)]] <- tibble(module=mod, p_AUPRC=p_aucpr, p_ROC=p_roc)
}

perm_all <- dplyr::bind_rows(perm_list)
pvals    <- dplyr::bind_rows(pval_list)

## =======================
## Plot: permutation boxes
## =======================
perm_long <- perm_all %>%
  dplyr::select(module, AUPRC, ROC) %>%
  tidyr::pivot_longer(c(AUPRC, ROC), names_to="Metric", values_to="Value")

obs_long <- report %>%
  dplyr::select(module, AUPRC, ROC) %>%
  tidyr::pivot_longer(c(AUPRC, ROC), names_to="Metric", values_to="Value") %>%
  dplyr::mutate(Type="Observed")

ggplot(perm_long, aes(x=factor(module), y=Value)) +
  geom_boxplot(fill="grey90", outlier.alpha=0.25, width=0.7) +
  geom_point(data=obs_long, aes(color=Metric), size=3, position=position_nudge(x=0.15)) +
  facet_wrap(~Metric, ncol=1, scales="free_y") +
  labs(x="Module", y="Score", title=sprintf("Permutation null (B=%d) with observed overlay", B)) +
  theme_classic(base_size=13) +
  theme(axis.title=element_text(face="bold"),
        axis.text=element_text(face="bold"))



## =======================
## Plot: permutation denisty
## =======================

library(dplyr)
library(ggplot2)
library(grid)
library(patchwork)   # install.packages("patchwork") if needed

## --- Subset to Modules 1,2,5,10 (keep your labeling) ---
perm_subset <- perm_long %>%
  filter(module %in% c(1,2,5,10), !is.na(Value)) %>%
  mutate(module_lab = factor(module, levels = c(1,2,5,10),
                             labels = c("Module 1","Module 2","Module 5","Module 10")))

obs_subset <- obs_long %>%
  filter(module %in% c(1,2,5,10), !is.na(Value)) %>%
  mutate(module_lab = factor(module, levels = c(1,2,5,10),
                             labels = c("Module 1","Module 2","Module 5","Module 10")))

# ---- Packages ----
library(dplyr)
library(ggplot2)
library(grid)
library(patchwork)  # to combine AUPRC + ROC
library(ggh4x)
load("RF+_RF-.RData")

# ---- Subset & labels (modules 1,2,5,10; adjust if needed) ----
perm_subset <- perm_long %>%
  filter(module %in% c(1,2,5,10), !is.na(Value)) %>%
  mutate(module_lab = factor(module, levels = c(1,2,5,10),
                             labels = c("Module 1","Module 2","Module 5","Module 10")))

obs_subset <- obs_long %>%
  filter(module %in% c(1,2,5,10), !is.na(Value)) %>%
  mutate(module_lab = factor(module, levels = c(1,2,5,10),
                             labels = c("Module 1","Module 2","Module 5","Module 10")))

# add once at top if not already
# install.packages("ggpattern")  # run once
library(ggpattern)
load("RF+_RF-.RData")

make_density <- function(metric, title_suffix = metric) {
  # observed per facet (this metric)
  obs_df <- obs_subset %>%
    dplyr::filter(Metric == metric) %>%
    dplyr::select(module_lab, obs = Value)
  
  # density curve per facet on a fixed grid
  base_dens <- perm_subset %>%
    dplyr::filter(Metric == metric) %>%
    dplyr::group_by(module_lab) %>%
    dplyr::reframe({
      d <- density(Value, adjust = 1.1, from = 0, to = 1.01, n = 512)
      tibble::tibble(x = d$x, y = pmax(d$y, 0))
    })
  
  # build a closed polygon for the area to the right of obs:
  # (trace along curve from obs->max, then return along baseline y=0 back to obs)
  poly_df <- base_dens %>%
    dplyr::left_join(obs_df, by = "module_lab") %>%
    dplyr::group_by(module_lab, obs) %>%
    dplyr::group_modify(function(df, key) {
      # keep right tail
      df_tail <- df[df$x >= key$obs, , drop = FALSE]
      if (nrow(df_tail) < 2) return(tibble::tibble())  # nothing to draw
      
      # ensure the polygon starts exactly at obs (linear interp for y at obs)
      if (df_tail$x[1] > key$obs) {
        # interpolate y at obs using the last point < obs and first >= obs
        left_idx  <- max(which(df$x < key$obs))
        right_idx <- min(which(df$x >= key$obs))
        if (is.finite(left_idx) && is.finite(right_idx)) {
          x0 <- df$x[left_idx]; y0 <- df$y[left_idx]
          x1 <- df$x[right_idx]; y1 <- df$y[right_idx]
          y_obs <- y0 + (y1 - y0) * ((key$obs - x0) / (x1 - x0))
          df_tail <- rbind(
            tibble::tibble(x = key$obs, y = max(y_obs, 0)),
            df_tail
          )
        }
      }
      
      # close the polygon along the baseline y=0 back to obs
      tibble::tibble(
        x = c(df_tail$x, rev(df_tail$x), df_tail$x[1]),
        y = c(df_tail$y, rep(0, nrow(df_tail)), df_tail$y[1] * 0)
      )
    }) %>%
    dplyr::ungroup()
  
  # p-values for labeling
  pval_df <- perm_subset %>%
    dplyr::filter(Metric == metric) %>%
    dplyr::left_join(obs_df, by = "module_lab") %>%
    dplyr::group_by(module_lab, obs) %>%
    dplyr::summarise(
      B_here = dplyr::n(),
      p = (sum(Value >= obs) + 1) / (B_here + 1),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      xlab = pmax(obs - 0.6 * obs, 0.0),
      lab  = sprintf("p = %.3g", p)
    )
  
  ggplot() +
    # --- hashed polygon under the curve for x >= obs ---
    ggpattern::geom_polygon_pattern(
      data = poly_df,
      aes(x = x, y = y, group = module_lab),
      inherit.aes = FALSE,
      pattern = "stripe",
      pattern_fill = "grey30",
      pattern_colour = NA,
      pattern_angle = 45,
      pattern_density = 0.5,
      pattern_spacing = 0.03,
      fill = NA, color = NA, alpha = 0
    ) +
    # full permutation density for context
    geom_density(
      data = dplyr::filter(perm_subset, Metric == metric),
      aes(x = Value, fill = "Permutation"),
      color = "black", alpha = 0.35, linewidth = 0.6, adjust = 1.1
    ) +
    # observed vertical line
    geom_vline(
      data = dplyr::filter(obs_subset, Metric == metric),
      aes(xintercept = Value, color = "Observed", linetype = "Observed"),
      linewidth = 0.8, show.legend = TRUE
    ) +
    # p-value label
    geom_text(
      data = pval_df,
      aes(x = xlab, y = Inf, label = lab),
      inherit.aes = FALSE,
      vjust = 1.1, fontface = "bold", size = 3.6
    ) +
    ggh4x::facet_wrap2(~ module_lab, scales = "free_y", axes = "all", strip.position = "top") +
    scale_x_continuous(limits = c(0, 1.01), breaks = seq(0, 1, 0.1),
                       expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Score", y = "Density",
         title = sprintf("Permutation density (B=%d) for %s (Modules 1,2,5 & 10)", B, title_suffix)) +
    scale_fill_manual(name = NULL, values = c("Permutation" = "green")) +
    scale_color_manual(name = NULL, values = c("Observed" = "red")) +
    scale_linetype_manual(name = NULL, values = c("Observed" = "dashed")) +
    guides(
      fill  = guide_legend(override.aes = list(color = NA)),
      color = guide_legend(override.aes = list(linetype = "dashed")),
      linetype = "none"
    ) +
    theme_classic(base_size = 11) +
    theme(
      panel.border   = element_rect(color = "black", linewidth = 0.6, fill = NA),
      strip.background = element_blank(),
      strip.text     = element_text(face = "bold"),
      axis.title     = element_text(face = "bold"),
      axis.text      = element_text(color = "black"),
      axis.ticks     = element_line(color = "black", linewidth = 0.4),
      axis.text.x    = element_text(size = 12, color = "black"),
      axis.text.y    = element_text(size = 13, color = "black"),
      axis.ticks.length = unit(3, "pt"),
      panel.spacing  = unit(10, "pt"),
      plot.title     = element_text(face = "bold", hjust = 0, margin = margin(b = 6)),
      legend.position = "top",
      legend.text    = element_text(size = 10, face = "bold")
    )
}




p_auprc <- make_density("AUPRC", "AUPRC")
p_roc   <- make_density("ROC",   "ROC")

ggsave(p_roc,file="/Users/javad/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/Figures/Final_presubmission_draft/Revision/Revision_Panels/RF+vsRF-_ROC.pdf",width = 10.1,height = 9.11,units = "in")
ggsave(p_auprc,file="/Users/javad/Library/CloudStorage/OneDrive-UniversityofPittsburgh/RA/Figures/Final_presubmission_draft/Revision/Revision_Panels/RF+vsRF-_AUPRC.pdf",width = 10.1,height = 9.11,units = "in")





