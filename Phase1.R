# ══════════════════════════════════════════════════════
# PHASE 1: Descriptive Statistics
# ══════════════════════════════════════════════════════

rm(list=ls())
library(zoo)
library(ggplot2)
library(moments)
library(MASS)
library(fBasics)
source("LoadFundData.R")

# ── Plots ─────────────────────────────────────────────
par(mfrow=c(2,1), cex=0.7, bty="l")
plot(P, main="S&P 500 Price Level",        xlab="", ylab="")
plot(r, main="S&P 500 Weekly Log Returns", xlab="", ylab="")

# ── Core variables ────────────────────────────────────
R   <- coredata(r)
N   <- length(R)
mu  <- mean(R)
sig <- sd(R)

# ── Moments (manual) ──────────────────────────────────
R0  <- R - mu
M2  <- sum(R0^2) / N
M3  <- sum(R0^3) / N
M4  <- sum(R0^4) / N
S        <- M3 / (sqrt(M2)^3)
K        <- M4 / (M2^2)        # total kurtosis
K_excess <- K - 3               # excess kurtosis

cat("=== Descriptive Statistics ===\n")
cat("Mean:          ", round(mu,       6), "\n")
cat("Std Dev:       ", round(sig,      6), "\n")
cat("Skewness:      ", round(S,        4), "\n")
cat("Kurtosis:      ", round(K,        4), "\n")
cat("Excess Kurtosis:", round(K_excess, 4), "\n")

# ── Annualized ────────────────────────────────────────
Nyear <- 365.25 / mean(as.numeric(diff(index(r))))
cat("\n=== Annualized ===\n")
cat("Annualized Mean:   ", round(mu  * Nyear,       4), "\n")
cat("Annualized Std Dev:", round(sig * sqrt(Nyear), 4), "\n")

# ── Normality Tests ───────────────────────────────────
cat("\n=== Normality Tests ===\n")
agostino.test(R)
anscombe.test(R)
jarque.test(R)

# ── Density comparison plot ───────────────────────────
x <- seq(-5, 5, 0.01)
df_plot <- data.frame(
  x    = rep(x, 4),
  y    = c(dnorm(x),
           dt(x, df=10) * sqrt(10/8),
           dt(x, df=5)  * sqrt(5/3),
           dt(x, df=3)  * sqrt(3/1)),
  dist = rep(c("Normal", "t (v=10)", "t (v=5)", "t (v=3)"), each=length(x))
)
ggplot(df_plot, aes(x=x, y=y, colour=dist)) +
  geom_line(linewidth=0.8) +
  labs(title="t-Student vs Normal distribution", y="Density", x="") +
  theme_bw() +
  scale_colour_manual("", values=c("black","red","blue","green"))

# ── QQ Plot — Normal ──────────────────────────────────
R0   <- (R - mu) / sig
q    <- seq(0.001, 0.999, 0.001)
Qemp <- quantile(R0, q)
Qteo <- qnorm(q)
lim0 <- c(-5, 5)

par(mfrow=c(1,1), cex=0.7, bty="l")
plot(Qteo, Qemp,
     main="QQ Plot (Normal)",
     col ="red",
     xlim=lim0, ylim=lim0,
     xlab="Theoretical Quantiles",
     ylab="Empirical Quantiles")
abline(a=0, b=1, lwd=2)

# ── Degrees of freedom — method of moments ────────────
v0 <- 4 + 6 / K_excess
cat("\nDegrees of freedom (moments):", round(v0, 3), "\n")

# ── Fit t-Student via MLE ─────────────────────────────
R0 <- as.numeric(scale(R))
d0 <- fitdistr(R0, "normal")
d1 <- fitdistr(R0, "t")
v  <- d1$estimate[["df"]]

cat("Degrees of freedom (MLE):    ", round(v, 3), "\n")
cat("Gain per obs (t vs normal):  ",
    round(100 * (d1$loglik - d0$loglik) / length(R0), 4), "\n")

# ── Histogram with fitted distributions (dùng v từ MLE) ──
R0_std <- (R - mu) / sig        # standardised returns
bwdth  <- 0.1
ggplot(data.frame(R0_std), aes(x=R0_std)) +
  theme_bw() +
  geom_histogram(binwidth=bwdth, colour="white",
                 fill="steelblue", linewidth=0.1) +
  stat_function(fun=function(x) dnorm(x) * N * bwdth,
                color="red",   linewidth=1) +
  stat_function(fun=function(x) dt(x, df=v) * sqrt((v-2)/v) * N * bwdth,
                color="black", linewidth=1) +
  labs(title=paste0("Empirical vs Normal vs t-Student (v=", round(v,2), ")"),
       x="Standardized Returns", y="Count")

# ── QQ Plot — t-Student (rescaled) ───────────────────
Qemp <- quantile(R0, q)
Qteo <- qt(q, df=v) * sqrt((v-2)/v)

plot(Qteo, Qemp,
     main="QQ Plot (t-Student)",
     col ="red",
     xlim=lim0, ylim=lim0,
     xlab="Theoretical Quantiles",
     ylab="Empirical Quantiles")
abline(a=0, b=1, lwd=2)