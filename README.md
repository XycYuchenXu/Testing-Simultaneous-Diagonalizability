# Testing Simultaneous Diagonalizability

This is the reproducible repository for the paper [*_Testing Simultaneous Diagonalizability_*](https://doi.org/10.1080/01621459.2023.2202435). It mainly relies on the `R` package [`eigTest`](https://github.com/XycYuchenXu/eigTest).

## Setup
For setup, first install some dependency packages (if not installed already) for reproduction with:
```
install.packages(c('devtools', 'foreach', 'parallel', 'doSNOW', 'tidyverse', 'tikzDevice', 'ggpattern', 'plotrix', 'reshape2', 'latex2exp', 'WDI', 'OECD', 'x12', 'pdfetch', 'jsonlite', 'MTS'))
```

The `R` package [`eigTest`](https://github.com/XycYuchenXu/eigTest) is available at Github. Brief usage descriptions are available there.

To install, run the following code in `R`:
```
devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = T)
```

For some implementation examples, please refer to the vignette with `browseVignettes('eigTest')`.

## Code
The folder [code](code) includes the scripts for different simulation studies and applications in the paper.

In particular, the simulation studies include examples for commutator-based test, LLR test, multi-sample eigenvector test and partial test. The exploratory script for two-sample tests under relatively high-dimensional setting is also provided. Applications include VAR examples on macroeconomic data of 8 countries, and Markovian model on Hudson river discharge data.

## Data
The folder [data](data) includes the `R` objects of:
1. the array of estimated VAR coefficient matrices for 8 countries `countryCoeff.RData`;
2. the estimated covariance of VAR coefficient matrices for 8 countries `countryCovar.RData`;
3. the weekly streamflow discharge data of Hudson river `hudsonWeekly.RData`;
4. the daily streamflow discharge data of Hudson river `hudsonDaily.RData`.
5. the replicates of raw p-values from simulation studies for test power analysis and histogram visualizations.

The first 4 objects are also available from the package [`eigTest`](https://github.com/XycYuchenXu/eigTest).

The only `csv` file [`M2_KOR_monthly_SA.csv`](data/M2_KOR_monthly_SA.csv) in the folder has the Korean money supply M2 data (source: Bank of Korea), which cannot be directly URL-accessed by functions in `R` from the website.

## Output
The folder [output](output) includes the replicates of raw p-values from simulation studies in folder [code](code), and the intermediate output of macroeconomic multivariate time series data `countryMacro.RData` accessed and processed by [`VAR_application.R`](code/VAR_application.R). The sub-folder [Plots](output/Plots) includes generated plots from both simulations and applications in [code](code), with either [`tikz`](output/Plots/tikz) or [`png`](output/Plots/png) format.
