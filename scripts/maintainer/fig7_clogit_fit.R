#!/usr/bin/env Rscript
# Conditional-logit fits behind Figure 7b-c and Supplementary Table SI_event_time_or_yearly (data/derived/fig7_prospective/or_clogit.tsv).
# Documentation of the model; not run by run_pipeline.sh. Its input, one row per (drug, decision year, candidate disease) with the
# outcome and the within-year standardized scores, is built from the licensed trial records and is not shipped.
# Strata = drug and decision year; candidates = diseases held out of benchmark tuning that are not yet trialed with the drug;
# outcome 1 = the drug's new trial diseases of that year. survival::clogit(method = "efron"), robust SEs clustered by drug.
# S1_* = one random case per stratum (single event, so breslow = exact), robust SEs.
# Usage: Rscript fig7_clogit_fit.R <data_dir> <out_dir> <years comma-sep> <tag> <mode yearly|pooled> <models comma-sep> [NNEG=0] [seed=1]
#   NNEG=0: full candidate set (used for all reported fits); NNEG>0: all positives + NNEG random negatives per stratum.
#   mode=pooled fits one model on all listed years together (set name "pooled_<min>_<max>"), strata still per drug and year.
# Reported fits: yearly 2000-2021 and pooled 2000-2021, models M0, M2, PERM_M2 and M2_<reference score>, sample "full", method "efron".
suppressPackageStartupMessages({library(survival); library(splines); library(data.table)})
a <- commandArgs(TRUE); DD <- a[1]; OUT <- a[2]; YRS <- as.integer(strsplit(a[3], ",")[[1]]); TAGARG <- a[4]; MODE <- a[5]
MLIST <- strsplit(a[6], ",")[[1]]; NNEG <- if (length(a) >= 7) as.integer(a[7]) else 0L; SEED <- if (length(a) >= 8) as.integer(a[8]) else 1L
FULLCHK <- FALSE; PSET <- sprintf("pooled_%d_%d", min(YRS), max(YRS))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); T0 <- Sys.time()
TAG <- TAGARG; F_OR <- file.path(OUT, sprintf("or_%s.tsv", TAG)); F_AU <- file.path(OUT, sprintf("auc_%s.tsv", TAG))
F_CN <- file.path(OUT, sprintf("counts_%s.tsv", TAG)); for (f in c(F_OR, F_AU, F_CN)) if (file.exists(f)) file.remove(f)
lg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf("[+%.0fs]", as.numeric(difftime(Sys.time(), T0, units = "secs"))), ..., "\n"); flush.console() }
COMP <- c("gene_jaccard", "adamic_adar_w", "gobp_simgic", "hpo_simgic", "strength")
rd <- function(y) fread(file.path(DD, sprintf("rows_y%d.tsv", y)),
                        select = c("year", "drug", "node", "y", "pop", paste0("z_", c("disatlas_ens5", COMP))))
FULL <- rbindlist(lapply(YRS, rd)); setnames(FULL, "z_disatlas_ens5", "z_dis")
FULL[, st := paste(year, drug)]; FULL[, lp := log1p(pop)]; FULL[, pop0 := as.integer(pop == 0)]
lg("full rows", nrow(FULL), "strata", uniqueN(FULL$st))
set.seed(SEED)
if (NNEG > 0) {
  FULL[, u := runif(.N)]
  FULL[, rk := frank(ifelse(y == 1, -1, u), ties.method = "first"), by = st]    # positives first, then random negatives
  FULL[, npos := sum(y), by = st]
  ALL <- FULL[y == 1 | rk <= npos + NNEG]
} else ALL <- FULL                                                                # full candidate set; no copy (memory)
ALL[, z_dis_perm := sample(z_dis), by = st]                                       # within-stratum permutation (sampled stratum)
ALL[, u1 := runif(.N)]; ALL[, keep1 := y == 0 | (y == 1 & u1 == min(u1[y == 1])), by = st]   # one random case per stratum
lg("sampled rows", nrow(ALL), "cases", sum(ALL$y))
for (s in if (MODE == "pooled") "pooled" else as.character(YRS)) {
  d <- if (s == "pooled") ALL else ALL[year == as.integer(s)]; k <- d[, .(np = sum(y), nn = sum(1 - y)), by = st]
  dF <- if (s == "pooled") FULL else FULL[year == as.integer(s)]
  fwrite(data.table(seed = SEED, set = if (s == "pooled") PSET else s, n_strata = nrow(k), n_rows_sampled = nrow(d), n_rows_full = nrow(dF), n_cases = sum(d$y),
                    n_drugs = uniqueN(d$drug), neg_per_stratum_min = min(k$nn), neg_per_stratum_median = median(k$nn),
                    strata_lt_NNEG_neg = sum(k$nn < NNEG), cases_per_stratum_median = median(k$np), cases_per_stratum_max = max(k$np),
                    strata_multi_case = sum(k$np > 1)), F_CN, sep = "\t", append = file.exists(F_CN))
}
# flexible popularity: natural cubic spline of log1p(pop) with 4 df + indicator pop==0. Most candidates have pop==0, so the
# default quantile knots of ns(df=4) all collapse onto 0; interior knots are therefore placed at the 25/50/75th percentiles of
# log1p(pop) among pop>0 rows of the fitted set (3 interior knots + 2 boundary knots = 4 df), the zeros being absorbed by pop0.
FLEX <- "NSPOP + pop0"
nsrhs <- function(rhs, d) { q <- quantile(d$lp[d$pop > 0], c(.25, .5, .75)); q <- unique(round(q, 6))
  gsub("NSPOP", sprintf("ns(lp, knots = c(%s), Boundary.knots = c(%s))", paste(q, collapse = ", "),
                        paste(round(range(d$lp), 6), collapse = ", ")), rhs, fixed = TRUE) }
MOD <- list(M0 = "z_dis", M1 = "z_dis + lp", M2 = paste("z_dis +", FLEX), P_log = "lp", P_flex = FLEX,
            PERM_M0 = "z_dis_perm", PERM_M2 = paste("z_dis_perm +", FLEX), H2H = paste("z_dis + z_gene_jaccard +", FLEX))
for (cm in COMP) { MOD[[paste0("M2_", cm)]] <- paste0("z_", cm, " + ", FLEX); MOD[[paste0("M0_", cm)]] <- paste0("z_", cm) }
auc_within <- function(d, lp) {
  x <- data.table(st = d$st, y = d$y, lp = lp)[, r := frank(lp, ties.method = "average"), by = st]
  s <- x[, .(np = sum(y), nn = sum(1 - y), rs = sum(r[y == 1])), by = st][np > 0 & nn > 0]
  s[, auc := (rs - np * (np + 1) / 2) / (np * nn)]
  data.table(auc_pairw = s[, sum(auc * np * nn) / sum(np * nn)], auc_stratmean = mean(s$auc))
}
fit1 <- function(setname, mname, rhs, method = "efron", sens = FALSE, full = FALSE) {
  src <- if (full) FULL else ALL
  d <- if (setname == "pooled") src else src[year == as.integer(setname)]
  if (sens) d <- d[keep1 == TRUE]
  rhs <- nsrhs(rhs, d); t0 <- Sys.time()
  # exact multi-case likelihood: survival gives no dfbeta for method="exact", so it is fitted without cluster (naive SEs)
  f <- tryCatch(if (method == "exact") clogit(as.formula(paste("y ~", rhs, "+ strata(st)")), data = d, method = "exact",
                                              control = coxph.control(iter.max = 50))
                else clogit(as.formula(paste("y ~", rhs, "+ strata(st)")), data = d, method = method, cluster = drug,
                            control = coxph.control(iter.max = 50)), error = function(e) e)
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  base <- data.table(seed = SEED, sample = if (full || NNEG == 0) "full" else sprintf("neg%d", NNEG), set = if (setname == "pooled") PSET else setname, model = mname, method = method,
                     n_strata = uniqueN(d$st), n_rows = nrow(d), n_cases = sum(d$y), n_drugs = uniqueN(d$drug), se_type = if (method == "exact") "naive" else "robust_cluster_drug", rhs = rhs, sec = round(sec, 1))
  if (inherits(f, "error")) { lg("ERROR", setname, mname, conditionMessage(f))
    fwrite(cbind(base, term = NA, error = conditionMessage(f)), F_OR, sep = "\t", append = file.exists(F_OR)); return(invisible()) }
  b <- coef(f); V <- vcov(f); Vn <- f$naive.var; if (is.null(Vn)) Vn <- V; se <- sqrt(diag(V))
  or <- cbind(base, data.table(term = names(b), coef = b, se_robust = se, se_naive = sqrt(diag(Vn)), OR = exp(b),
                               lo95 = exp(b - 1.959964 * se), hi95 = exp(b + 1.959964 * se), p_robust = 2 * pnorm(-abs(b / se)),
                               loglik0 = f$loglik[1], loglik = f$loglik[2], converged = f$iter < 50, error = NA_character_))
  fwrite(or, F_OR, sep = "\t", append = file.exists(F_OR))
  fwrite(cbind(base, auc_within(d, predict(f, type = "lp"))), F_AU, sep = "\t", append = file.exists(F_AU))
  zt <- grep("^z_", names(b), value = TRUE)
  lg(setname, mname, sprintf("%.1fs", sec), if (length(zt)) paste(sprintf("%s OR=%.3f [%.3f,%.3f]", zt, exp(b[zt]), exp(b[zt] - 1.96 * se[zt]), exp(b[zt] + 1.96 * se[zt])), collapse = "; ") else "")
}
Y <- as.character(YRS)
run1 <- function(s, m) if (startsWith(m, "S1_")) fit1(s, m, MOD[[sub("^S1_", "", m)]], method = "breslow", sens = TRUE) else fit1(s, m, MOD[[m]])
stopifnot(all(sub("^S1_", "", MLIST) %in% names(MOD)))
if (MODE == "pooled") for (m in MLIST) run1("pooled", m) else for (s in Y) for (m in MLIST) run1(s, m)
lg("done")
