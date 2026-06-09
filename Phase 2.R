##########################################################
# Phase 2 — Unconditional VaR & ES + Backtesting        #
##########################################################

rm(list=ls())
require(zoo)
require(rugarch)
require(MASS)

source("LoadFundData.R")

R   <- coredata(r)
N   <- length(R)
mu  <- mean(R)
sig <- sd(R)

w_length <- 250
rr       <- tail(r, N - w_length)

cat("Observations:        ", N,                    "\n")
cat("Out-of-sample obs:   ", length(rr),           "\n")
cat("Mean:                ", round(mu,  6),        "\n")
cat("Std Dev:             ", round(sig, 6),        "\n")

##########################################################
# STEP 1 — Compute VaR & ES                             #
##########################################################

alpha_95 <- 0.05
alpha_99 <- 0.01

# ── Normal ────────────────────────────────────────────
VaR_N_95 <- mu + sig * qnorm(alpha_95)
VaR_N_99 <- mu + sig * qnorm(alpha_99)
ES_N_95  <- mu - sig * dnorm(qnorm(alpha_95)) / alpha_95
ES_N_99  <- mu - sig * dnorm(qnorm(alpha_99)) / alpha_99

cat("\n=== Normal ===\n")
cat("VaR 95%:", round(VaR_N_95, 6), "\n")
cat("VaR 99%:", round(VaR_N_99, 6), "\n")
cat("ES  95%:", round(ES_N_95,  6), "\n")
cat("ES  99%:", round(ES_N_99,  6), "\n")

# ── t-Student ─────────────────────────────────────────
R0 <- as.numeric(scale(R))
d0 <- fitdistr(R0, "normal")
d1 <- fitdistr(R0, "t")
v  <- d1$estimate[["df"]]
cat("\nDegrees of freedom:", round(v, 3), "\n")

VaR_t_95 <- mu + sig * qt(alpha_95, df=v) * sqrt((v-2)/v)
VaR_t_99 <- mu + sig * qt(alpha_99, df=v) * sqrt((v-2)/v)
ES_t_95  <- mu - sig * (dt(qt(alpha_95,df=v),df=v)/alpha_95) *
  ((v+qt(alpha_95,df=v)^2)/(v-1)) * sqrt((v-2)/v)
ES_t_99  <- mu - sig * (dt(qt(alpha_99,df=v),df=v)/alpha_99) *
  ((v+qt(alpha_99,df=v)^2)/(v-1)) * sqrt((v-2)/v)

cat("\n=== t-Student ===\n")
cat("Degrees of freedom:", round(v,       3), "\n")
cat("VaR 95%:",           round(VaR_t_95, 6), "\n")
cat("VaR 99%:",           round(VaR_t_99, 6), "\n")
cat("ES  95%:",           round(ES_t_95,  6), "\n")
cat("ES  99%:",           round(ES_t_99,  6), "\n")

# ── Historical Simulation ─────────────────────────────
VaR_HS_95 <- quantile(R, probs=alpha_95)
VaR_HS_99 <- quantile(R, probs=alpha_99)
ES_HS_95  <- mean(R[R <= VaR_HS_95])
ES_HS_99  <- mean(R[R <= VaR_HS_99])

cat("\n=== Historical Simulation ===\n")
cat("VaR 95%:", round(VaR_HS_95, 6), "\n")
cat("VaR 99%:", round(VaR_HS_99, 6), "\n")
cat("ES  95%:", round(ES_HS_95,  6), "\n")
cat("ES  99%:", round(ES_HS_99,  6), "\n")

##########################################################
# STEP 2 — Backtesting Function                         #
##########################################################

backtest <- function(label, p, rr, var_bt, es_bt) {
  cat("\n=== Backtesting:", label, (1-p)*100, "% ===\n")
  
  temp  <- VaRTest(alpha=p,
                   actual=coredata(rr),
                   VaR=coredata(var_bt))
  
  temp1 <- VaRDurTest(alpha=p,
                      actual=coredata(rr),
                      VaR=coredata(var_bt),
                      conf.level=0.95)
  
  # Ch2 independence = LRcc - LRuc
  LRind_stat <- temp$cc.LRstat - temp$uc.LRstat
  LRind_p    <- pchisq(LRind_stat, df=1, lower.tail=FALSE)
  
  cat("Expected violations:", round(length(rr)*p), "\n")
  cat("Actual violations:  ", temp$actual.exceed,  "\n")
  cat("Kupiec  p-value:", round(temp$uc.LRp,  4), "\n")
  cat("Ch2     p-value:", round(LRind_p,      4), "\n")
  cat("Ch-P    p-value:", round(temp1$LRp,    4), "\n")
  
  es_res <- ESTest(alpha=p,
                   actual=coredata(rr),
                   ES=coredata(es_bt),
                   VaR=coredata(var_bt))
  cat("ES      p-value:", round(es_res$p.value, 4), "\n")
  cat("ES Decision:",     es_res$Decision,          "\n")
  
  VaRplot(alpha=p, actual=rr, VaR=var_bt)
  title(main=paste0(label, ": VaR violations (", (1-p)*100, "%)"))
}

##########################################################
# STEP 3 — Run Backtesting at 95%                       #
##########################################################

p <- 0.05

varN_bt  <- zoo(rep(VaR_N_95,  length(rr)), order.by=index(rr))
esN_bt   <- zoo(rep(ES_N_95,   length(rr)), order.by=index(rr))
varT_bt  <- zoo(rep(VaR_t_95,  length(rr)), order.by=index(rr))
esT_bt   <- zoo(rep(ES_t_95,   length(rr)), order.by=index(rr))
varHS_bt <- zoo(rep(VaR_HS_95, length(rr)), order.by=index(rr))
esHS_bt  <- zoo(rep(ES_HS_95,  length(rr)), order.by=index(rr))

backtest("Normal",     p, rr, varN_bt,  esN_bt)
backtest("t-Student",  p, rr, varT_bt,  esT_bt)
backtest("Historical", p, rr, varHS_bt, esHS_bt)

##########################################################
# STEP 4 — Run Backtesting at 99%                       #
##########################################################

p <- 0.01

varN_bt  <- zoo(rep(VaR_N_99,  length(rr)), order.by=index(rr))
esN_bt   <- zoo(rep(ES_N_99,   length(rr)), order.by=index(rr))
varT_bt  <- zoo(rep(VaR_t_99,  length(rr)), order.by=index(rr))
esT_bt   <- zoo(rep(ES_t_99,   length(rr)), order.by=index(rr))
varHS_bt <- zoo(rep(VaR_HS_99, length(rr)), order.by=index(rr))
esHS_bt  <- zoo(rep(ES_HS_99,  length(rr)), order.by=index(rr))

backtest("Normal",     p, rr, varN_bt,  esN_bt)
backtest("t-Student",  p, rr, varT_bt,  esT_bt)
backtest("Historical", p, rr, varHS_bt, esHS_bt)

##########################################################
# Stress Test — Unconditional Models                    #
##########################################################

# ── Define periods ────────────────────────────────────
periods <- list(
  "Calm 2013-2019"  = c("2013-01-01", "2019-12-31"),
  "COVID 2020"      = c("2020-01-01", "2020-12-31"),
  "Rate Hike 2022"  = c("2022-01-01", "2022-12-31"),
  "Post 2022"       = c("2023-01-01", "2026-03-25"),
  "Full Sample"     = c("2012-04-01", "2026-03-25")
)

# ── Stress test function ──────────────────────────────
stress_period <- function(period_name, dates, rr, models_list, p) {
  cat("\n--- Period:", period_name, "---\n")
  start_d <- as.Date(dates[1])
  end_d   <- as.Date(dates[2])
  rr_sub  <- window(rr, start=start_d, end=end_d)
  
  if(length(rr_sub) < 5) {
    cat("Not enough observations\n")
    return()
  }
  
  cat(sprintf("%-25s %4s %6s %8s\n",
              "Model", "Obs", "Viol", "Rate(%)"))
  cat(strrep("-", 50), "\n")
  
  for(model_name in names(models_list)) {
    var_sub  <- window(models_list[[model_name]],
                       start=start_d, end=end_d)
    viol     <- sum(coredata(rr_sub) < coredata(var_sub),
                    na.rm=TRUE)
    total    <- sum(!is.na(coredata(rr_sub)))
    rate     <- round(viol / total * 100, 1)
    expected <- round(p * 100, 1)
    flag     <- ifelse(abs(rate - expected) > expected * 0.5,
                       "<<", "")
    cat(sprintf("%-25s %4d %6d %7.1f%% %s\n",
                model_name, total, viol, rate, flag))
  }
}

# ── Build model lists ─────────────────────────────────
varN_99_bt  <- zoo(rep(VaR_N_99,  length(rr)), order.by=index(rr))
varT_99_bt  <- zoo(rep(VaR_t_99,  length(rr)), order.by=index(rr))
varHS_99_bt <- zoo(rep(VaR_HS_99, length(rr)), order.by=index(rr))

models_95 <- list(
  "Normal 95%"     = varN_bt,
  "t-Student 95%"  = varT_bt,
  "Historical 95%" = varHS_bt
)

models_99 <- list(
  "Normal 99%"     = varN_99_bt,
  "t-Student 99%"  = varT_99_bt,
  "Historical 99%" = varHS_99_bt
)

# ── Run stress test ───────────────────────────────────
cat("\n========================================\n")
cat("STRESS TEST — 95% VaR\n")
cat("Expected violation rate: 5%\n")
cat("========================================\n")

for(period_name in names(periods)) {
  stress_period(period_name, periods[[period_name]],
                rr, models_95, p=0.05)
}

cat("\n========================================\n")
cat("STRESS TEST — 99% VaR\n")
cat("Expected violation rate: 1%\n")
cat("========================================\n")

for(period_name in names(periods)) {
  stress_period(period_name, periods[[period_name]],
                rr, models_99, p=0.01)
}