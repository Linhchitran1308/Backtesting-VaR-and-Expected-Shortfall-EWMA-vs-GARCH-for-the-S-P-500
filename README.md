# Backtesting VaR and Expected Shortfall: EWMA vs GARCH for the S&P 500

**Bachelor's Thesis — SGH Warsaw School of Economics, 2026**  
**Author:** Linh Chi Tran | **Supervisor:** Prof. dr hab. Michał Rubaszek

---

## Research Question

Which combination of volatility model (EWMA or GARCH) and distributional assumption (Normal or t-Student) produces the most statistically adequate Value at Risk (VaR) and Expected Shortfall (ES) estimates for the S&P 500, as evaluated by formal regulatory backtesting at the 95% and 99% confidence levels?

---

## Key Finding

> **GARCH(1,1) with t-Student innovations is the only specification to pass all four regulatory backtesting criteria simultaneously at the 99% confidence level** — making it the recommended model for Basel III internal-models capital reporting.

---

## Methodology

### Models
| Model | Type | Key Property |
|-------|------|-------------|
| Unconditional Normal | Static | Full-sample variance |
| Unconditional t-Student | Static | Full-sample, fat tails (ν = 9.17) |
| Historical Simulation | Static | Empirical quantile |
| EWMA (λ = 0.94) | Dynamic | Geometrically declining weights, no mean reversion |
| GARCH(1,1) Normal | Dynamic | MLE, reverts to equilibrium variance |
| GARCH(1,1) t-Student | Dynamic | MLE, fat-tailed innovations |

### Backtesting Framework
Four formal tests applied to all 14 model–distribution combinations. **PASS = all four fail to reject at 5% simultaneously.**

| Test | H₀ |
|------|----|
| Kupiec (1995) | Correct violation frequency |
| Christoffersen Ch2 (1998) | Violations are independent |
| Christoffersen-Pelletier Ch-P (2004) | Durations between violations are memoryless |
| McNeil-Frey ES (2000) | Tail loss magnitude is adequate |

### Data
- **Index:** S&P 500 weekly log-returns (`^spx` from stooq.pl)
- **Sample:** January 2007 – March 2026 (~1,000 observations)
- **Out-of-sample:** 752 weeks (rolling 250-week estimation window)

---

## Results Summary

| Model | Level | Kupiec | Ch2 | Ch-P | ES | Overall |
|-------|-------|--------|-----|------|----|---------|
| EWMA Normal | 95% | ✓ | ✓ | ✓ | ✗ | PARTIAL |
| EWMA Normal | 99% | ✗ | ✗ | ✗ | ✗ | FAIL |
| EWMA t-Student | 95% | ✓ | ✓ | ✓ | ✗ | PARTIAL |
| EWMA t-Student | 99% | ✗ | ✗ | ✗ | ✗ | FAIL |
| GARCH Normal | 95% | ✗ | ✗ | ✓ | ✓ | PARTIAL |
| GARCH Normal | 99% | ✗ | ✗ | ✓ | ✗ | PARTIAL |
| GARCH t-Student | 95% | ✗ | ✗ | ✓ | ✓ | PARTIAL |
| **GARCH t-Student** | **99%** | **✓** | **✓** | **✓** | **✓** | **PASS ⭐** |

### Three Cross-Cutting Patterns
1. **ES is the hardest test** — fails every model at 95% except Historical and GARCH t-Student
2. **GARCH always passes Ch-P** — reversion to equilibrium variance prevents violation clustering regardless of distribution
3. **Model selection is confidence-level specific** — EWMA better calibrated at 95%, GARCH t-Student is the only adequate specification at 99%

### Stress Test Highlights (95% VaR, target: 5%)
| Period | Unconditional | EWMA Normal | GARCH t-Student |
|--------|--------------|-------------|-----------------|
| Calm 2013–2019 | 0.3–1.1% ❌ | 4.9% ✓ | 6.8% ✓ |
| COVID-19 2020 | 5.8–7.7% | 7.7% | 9.6% |
| Rate Hike 2022 | 0.0–3.8% | 5.8% ✓ | 7.7% |
| Full sample | 0.7–1.5% ❌ | 4.4% ✓ | 7.3% |

---

## Repository Structure

```
├── LoadFundData.R               # Data loading (S&P 500 from stooq.pl)
├── Phase1_DescriptiveStats.R    # Descriptive statistics, QQ plots, MLE fit
├── Phase2_UnconditionalVaR.R    # Unconditional VaR/ES + backtesting
├── Phase3_VolatilityModels.R    # EWMA & GARCH rolling window + backtesting
└── output/
    └── LinhChiTran_MR_thesis.pdf
```

---

## How to Run

```r
# 1. Install required packages
install.packages(c("zoo", "rugarch", "MASS", "ggplot2",
                   "moments", "fBasics"))

# 2. Load data (requires internet connection to stooq.pl)
source("LoadFundData.R")

# 3. Run phases in order
source("Phase1_DescriptiveStats.R")
source("Phase2_UnconditionalVaR.R")
source("Phase3_VolatilityModels.R")  # ~10–15 min runtime (rolling window)
```

> **Note:** Phase 3 runs ~1,000 rolling window iterations and may take 10–15 minutes depending on hardware.

---

## Tech Stack

![R](https://img.shields.io/badge/R-4.x-276DC3?logo=r&logoColor=white)
![rugarch](https://img.shields.io/badge/rugarch-GARCH%20models-blue)
![ggplot2](https://img.shields.io/badge/ggplot2-visualisation-red)

| Package | Purpose |
|---------|---------|
| `rugarch` | GARCH/EWMA specification, fitting, forecasting, backtesting |
| `zoo` | Time series handling and rolling windows |
| `MASS` | t-Student MLE via `fitdistr` |
| `ggplot2` | QQ plots, histograms, density comparisons |
| `moments` / `fBasics` | Skewness, kurtosis, normality tests |

---

## Citation

> Tran, L.C. (2026). *Backtesting VaR and Expected Shortfall: A Comparison of EWMA and GARCH Models for the S&P 500 Index.* Bachelor's thesis, SGH Warsaw School of Economics.
