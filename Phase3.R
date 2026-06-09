##########################################################
# Phase 3 — EWMA & GARCH Rolling VaR/ES                #
# Backtesting + Stress Test — Dynamic Models            #
##########################################################

rm(list=ls())
require(zoo)
require(rugarch)

source("LoadFundData.R")

R        <- coredata(r)
N        <- length(R)
w_length <- 250
rr       <- tail(r, N - w_length)
lambda   <- 0.94

##########################################################
# Helper Functions                                      #
##########################################################

es_normal <- function(sigma, p) {
  -sigma * dnorm(qnorm(p)) / p
}

es_tstudent <- function(sigma, shape, p) {
  q <- qt(p, df=shape)
  -sigma * sqrt((shape-2)/shape) *
    (dt(q, df=shape) / p) *
    ((shape + q^2) / (shape - 1))
}

print_backtest <- function(label, p, rr, var_bt, es_bt) {
  temp  <- VaRTest(alpha=p,
                   actual=coredata(rr),
                   VaR=coredata(var_bt))
  temp1 <- VaRDurTest(alpha=p,
                      actual=coredata(rr),
                      VaR=coredata(var_bt),
                      conf.level=0.95)
  LRind_stat <- temp$cc.LRstat - temp$uc.LRstat
  LRind_p    <- pchisq(LRind_stat, df=1, lower.tail=FALSE)
  es_res     <- ESTest(alpha=p,
                       actual=coredata(rr),
                       ES=coredata(es_bt),
                       VaR=coredata(var_bt))
  
  cat("\n===", label, (1-p)*100, "% ===\n")
  cat("Expected violations:", round(length(rr)*p), "\n")
  cat("Actual violations:  ", temp$actual.exceed,  "\n")
  cat("Kupiec  p-value:", round(temp$uc.LRp,    4), "\n")
  cat("Ch2     p-value:", round(LRind_p,        4), "\n")
  cat("Ch-P    p-value:", round(temp1$LRp,      4), "\n")
  cat("ES      p-value:", round(es_res$p.value, 4), "\n")
  cat("ES Decision:    ", es_res$Decision,          "\n")
  
  VaRplot(alpha=p, actual=rr, VaR=var_bt)
  title(main=paste0(label, ": VaR violations (", (1-p)*100, "%)"))
}

stress_period <- function(period_name, dates, rr,
                          models_list, p) {
  cat("\n--- Period:", period_name, "---\n")
  start_d <- as.Date(dates[1])
  end_d   <- as.Date(dates[2])
  rr_sub  <- window(rr, start=start_d, end=end_d)
  
  if (length(rr_sub) < 5) {
    cat("Not enough observations\n")
    return()
  }
  
  cat(sprintf("%-25s %4s %6s %8s\n",
              "Model", "Obs", "Viol", "Rate(%)"))
  cat(strrep("-", 50), "\n")
  
  for (model_name in names(models_list)) {
    var_sub <- window(models_list[[model_name]],
                      start=start_d, end=end_d)
    valid   <- !is.na(coredata(var_sub))
    rr_v    <- coredata(rr_sub)[valid]
    var_v   <- coredata(var_sub)[valid]
    viol    <- sum(rr_v < var_v)
    total   <- length(rr_v)
    if (total < 1) { cat(model_name, ": no data\n"); next }
    rate    <- round(viol / total * 100, 1)
    flag    <- ifelse(abs(rate - p*100) > p*100*0.5, "<<", "")
    cat(sprintf("%-25s %4d %6d %7.1f%% %s\n",
                model_name, total, viol, rate, flag))
  }
}

periods <- list(
  "Calm 2013-2019"  = c("2013-01-01", "2019-12-31"),
  "COVID 2020"      = c("2020-01-01", "2020-12-31"),
  "Rate Hike 2022"  = c("2022-01-01", "2022-12-31"),
  "Post 2022"       = c("2023-01-01", "2026-03-25"),
  "Full Sample"     = c("2012-04-01", "2026-03-25")
)

##########################################################
# SPECS                                                 #
##########################################################

shape_uncond <- 9.166

EWMAspec_N <- ugarchspec(
  mean.model         = list(armaOrder=c(0,0), include.mean=FALSE),
  variance.model     = list(model="iGARCH", garchOrder=c(1,1)),
  fixed.pars         = list(alpha1=1-lambda, omega=0),
  distribution.model = "norm")

EWMAspec_t <- ugarchspec(
  mean.model         = list(armaOrder=c(0,0), include.mean=FALSE),
  variance.model     = list(model="iGARCH", garchOrder=c(1,1)),
  fixed.pars         = list(alpha1=1-lambda, omega=0, shape=shape_uncond),
  distribution.model = "std")

GARCHspec_N <- ugarchspec(
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  variance.model     = list(model="sGARCH", garchOrder=c(1,1)),
  distribution.model = "norm")

GARCHspec_t <- ugarchspec(
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  variance.model     = list(model="sGARCH", garchOrder=c(1,1)),
  distribution.model = "std")

##########################################################
# EWMA Normal — rolling window                         #
##########################################################

cat("Running EWMA Normal...\n")
vares_EN <- rollapply(r, width=w_length,
                      function(w) {
                        frc   <- ugarchforecast(EWMAspec_N, data=w, n.ahead=1)
                        sigma <- as.numeric(sigma(frc))
                        var95 <- as.numeric(quantile(frc, 0.05))
                        var99 <- as.numeric(quantile(frc, 0.01))
                        es95  <- es_normal(sigma, 0.05)
                        es99  <- es_normal(sigma, 0.01)
                        return(c(var95, es95, var99, es99))
                      }, by=1, align="right")
vares_EN     <- stats::lag(vares_EN, -1)
varEWMA_N95  <- vares_EN[, 1];  esEWMA_N95  <- vares_EN[, 2]
varEWMA_N99  <- vares_EN[, 3];  esEWMA_N99  <- vares_EN[, 4]
cat("Done!\n")

# Backtesting
for (p in c(0.05, 0.01)) {
  var_bt <- if (p == 0.05) varEWMA_N95 else varEWMA_N99
  es_bt  <- if (p == 0.05) esEWMA_N95  else esEWMA_N99
  print_backtest("EWMA Normal", p, rr, var_bt, es_bt)
}

##########################################################
# EWMA t-Student — rolling window                      #
##########################################################

cat("Running EWMA t-Student...\n")
vares_Et <- rollapply(r, width=w_length,
                      function(w) {
                        frc   <- ugarchforecast(EWMAspec_t, data=w, n.ahead=1)
                        sigma <- as.numeric(sigma(frc))
                        var95 <- as.numeric(quantile(frc, 0.05))
                        var99 <- as.numeric(quantile(frc, 0.01))
                        es95  <- es_tstudent(sigma, shape_uncond, 0.05)
                        es99  <- es_tstudent(sigma, shape_uncond, 0.01)
                        return(c(var95, es95, var99, es99))
                      }, by=1, align="right")
vares_Et     <- stats::lag(vares_Et, -1)
varEWMA_t95  <- vares_Et[, 1];  esEWMA_t95  <- vares_Et[, 2]
varEWMA_t99  <- vares_Et[, 3];  esEWMA_t99  <- vares_Et[, 4]
cat("Done!\n")

for (p in c(0.05, 0.01)) {
  var_bt <- if (p == 0.05) varEWMA_t95 else varEWMA_t99
  es_bt  <- if (p == 0.05) esEWMA_t95  else esEWMA_t99
  print_backtest("EWMA t-Student", p, rr, var_bt, es_bt)
}

##########################################################
# GARCH Normal — rolling window                        #
##########################################################

cat("Running GARCH Normal... (~5 min)\n")
vares_GN <- rollapply(r, width=w_length,
                      function(w) {
                        fit   <- ugarchfit(data=w, spec=GARCHspec_N, solver="hybrid")
                        frc   <- ugarchforecast(fit, n.ahead=1)
                        sigma <- as.numeric(sigma(frc))
                        var95 <- as.numeric(quantile(frc, 0.05))
                        var99 <- as.numeric(quantile(frc, 0.01))
                        es95  <- es_normal(sigma, 0.05)
                        es99  <- es_normal(sigma, 0.01)
                        return(c(var95, es95, var99, es99))
                      }, by=1, align="right")
vares_GN     <- stats::lag(vares_GN, -1)
varGARCH_N95 <- vares_GN[, 1];  esGARCH_N95 <- vares_GN[, 2]
varGARCH_N99 <- vares_GN[, 3];  esGARCH_N99 <- vares_GN[, 4]
cat("Done!\n")

for (p in c(0.05, 0.01)) {
  var_bt <- if (p == 0.05) varGARCH_N95 else varGARCH_N99
  es_bt  <- if (p == 0.05) esGARCH_N95  else esGARCH_N99
  print_backtest("GARCH Normal", p, rr, var_bt, es_bt)
}

##########################################################
# GARCH t-Student — rolling window                     #
##########################################################

cat("Running GARCH t-Student... (~5 min)\n")
vares_Gt <- rollapply(r, width=w_length,
                      function(w) {
                        fit   <- ugarchfit(data=w, spec=GARCHspec_t, solver="hybrid")
                        frc   <- ugarchforecast(fit, n.ahead=1)
                        sigma <- as.numeric(sigma(frc))
                        shape <- coef(fit)["shape"]
                        var95 <- as.numeric(quantile(frc, 0.05))
                        var99 <- as.numeric(quantile(frc, 0.01))
                        es95  <- es_tstudent(sigma, shape, 0.05)
                        es99  <- es_tstudent(sigma, shape, 0.01)
                        return(c(var95, es95, var99, es99))
                      }, by=1, align="right")
vares_Gt     <- stats::lag(vares_Gt, -1)
varGARCH_t95 <- vares_Gt[, 1];  esGARCH_t95 <- vares_Gt[, 2]
varGARCH_t99 <- vares_Gt[, 3];  esGARCH_t99 <- vares_Gt[, 4]
cat("Done!\n")

for (p in c(0.05, 0.01)) {
  var_bt <- if (p == 0.05) varGARCH_t95 else varGARCH_t99
  es_bt  <- if (p == 0.05) esGARCH_t95  else esGARCH_t99
  print_backtest("GARCH t-Student", p, rr, var_bt, es_bt)
}

##########################################################
# GARCH Parameter Estimates — Table 3.3                #
##########################################################

fit_n  <- ugarchfit(data=r, spec=GARCHspec_N, solver="hybrid")
fit_t  <- ugarchfit(data=r, spec=GARCHspec_t, solver="hybrid")
coef_n <- coef(fit_n)
coef_t <- coef(fit_t)

cat("\n========================================\n")
cat("TABLE 3.3 — GARCH Parameter Estimates\n")
cat("========================================\n")
cat(sprintf("%-25s %12s %16s\n",
            "Parameter", "GARCH Normal", "GARCH t-Student"))
cat(strrep("-", 55), "\n")
cat(sprintf("%-25s %12.6f %16.6f\n", "omega",
            coef_n["omega"],  coef_t["omega"]))
cat(sprintf("%-25s %12.4f %16.4f\n", "alpha1 (ARCH)",
            coef_n["alpha1"], coef_t["alpha1"]))
cat(sprintf("%-25s %12.4f %16.4f\n", "beta1 (GARCH)",
            coef_n["beta1"],  coef_t["beta1"]))
cat(sprintf("%-25s %12.4f %16.4f\n", "alpha+beta",
            coef_n["alpha1"]+coef_n["beta1"],
            coef_t["alpha1"]+coef_t["beta1"]))
cat(sprintf("%-25s %12s %16.4f\n",   "nu (shape)",
            "—", coef_t["shape"]))
cat(sprintf("%-25s %12.2f %16.2f\n", "Half-life (weeks)",
            log(0.5)/log(coef_n["alpha1"]+coef_n["beta1"]),
            log(0.5)/log(coef_t["alpha1"]+coef_t["beta1"])))

##########################################################
# Stress Test — Dynamic Models                         #
##########################################################

models_95 <- list(
  "EWMA Normal"     = varEWMA_N95,
  "EWMA t-Student"  = varEWMA_t95,
  "GARCH Normal"    = varGARCH_N95,
  "GARCH t-Student" = varGARCH_t95
)

models_99 <- list(
  "EWMA Normal"     = varEWMA_N99,
  "EWMA t-Student"  = varEWMA_t99,
  "GARCH Normal"    = varGARCH_N99,
  "GARCH t-Student" = varGARCH_t99
)

cat("\n========================================\n")
cat("STRESS TEST — DYNAMIC MODELS — 95% VaR\n")
cat("Expected violation rate: 5%\n")
cat("========================================\n")
for (period_name in names(periods)) {
  stress_period(period_name, periods[[period_name]],
                rr, models_95, p=0.05)
}

cat("\n========================================\n")
cat("STRESS TEST — DYNAMIC MODELS — 99% VaR\n")
cat("Expected violation rate: 1%\n")
cat("========================================\n")
for (period_name in names(periods)) {
  stress_period(period_name, periods[[period_name]],
                rr, models_99, p=0.01)
}
Tóm tắt thay đổi so với Phase 3 + Phase 4 cũ:
  
  Gộp EWMA Normal vào — trước đây Phase 3 thiếu
Mỗi model chạy 1 lần cho cả 95% và 99% — tiết kiệm ~50% runtime
Stress test dynamic ở cuối file — không cần Phase 4 riêng
Bỏ Phase 4 hoàn toàn — repo giờ chỉ còn 3 phases sạch gọn

Repo structure cuối cùng:
  ├── README.md
├── LoadFundData.R
├── Phase1_DescriptiveStats.R
├── Phase2_UnconditionalVaR.R
├── Phase3_VolatilityModels.R     ← file này
└── output/
  └── LinhChiTran_MR_thesis.pdf