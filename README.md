# Testing Simultaneous Diagonalizability

This is the reproducible repository for the paper. It mainly relies on the `R` package [`eigTest`](https://github.com/XycYuchenXu/eigTest) developed by the authors.

## Setup
The `R` package [`eigTest`](https://github.com/XycYuchenXu/eigTest) is available at Github for testing simultaneous diagonalizability.

To install, run the following code in `R`:
`devtools::install_github('XycYuchenXu/eigTest', force = T, build_vignettes = T)`

## Usage
That package has the following functionalities:

1. Check whether the means of two random square matrices (`matrixA` and `matrixB`) with dimension `d`-by-`d` are simultaneously diagonalizable, either considering the commutator of the two matrices (use function `commutatorTest`), or using the log-likelihood ratio test framework (`projTest`), given asymptotic limiting covariance matrices (`covMatA` and `covMatB`) and a convergence rate (`cn`). For both functions, one needs to input an array of two matrices (`A` such that `A[1,,] = matrixA` and `A[2,,] = matrixB`), an array of two limiting covariance matrices (`cov.arr` such that `cov.arr[1,,] = covMatA` and `cov.arr[2,,] = covMatB`), and the convergence rate (`cn`).

2. Given an array of square random matrices (`A`), estimate a common eigenvector matrix `V` (by supplying the array of matrices `A` for function `JDTE`).

3. Check whether a set of matrices (by testing on the estimator array `A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) are simultaneously diagonalizable. The test is conducted by checking whether the eigenvector matrix estimated above (`V`, the output from function `JDTE`) is indeed the common eigenvector matrix. Use function `eigTest` with input of two arrays (`A` and `cov.arr`) and the convergence rate `cn`. `V` is optional as `JDTE` will be called by default internally.

4. Given an array of square random matrices (`A` with dimension `p`-by-`d`-by-`d`), estimate `k`-dimensional (`k` ≤ `d`) common invariant sub-space `Qk` (by supplying the array of matrices `A` and parameter `k` for function `expmPartSchur`). The output `Q` from `expmPartSchur` satisfies `Q[,1:k] = Qk`.

5. Check whether a set of matrices (by testing on the estimator array `A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) share `k` of the all `d` eigenvectors (`k` ≤ `d`). The test is conducted by checking whether the estimated `k`-dimensional invariant sub-space `Qk` from `expmPartSchur` and the estimated common eigenvectors `V` for `B`, where `B[i,,] = Qk' A[i,,] Qk`, satisfy structure properties. Use function `partialTest` with input of two arrays (`A` and `cov.arr`), the convergence rate `cn`, and number of common eigenvectors `k`. `Qk` and `V` are optional as `expmPartSchur` and `JDTE` will be called by default internally.

6. Generate Gaussian samples for simulations. The function `generateMeans` is used to generate `p` mean matrices with dimension `d`-by-`d` and noise levels `snr` on common eigenvector `V`, the output is an array of dimension `p`-by-`q`-by-`d`-by-`d` where `q = length(snr) + 1`. Function `simuSamples` takes the output from `generateMeans` as input for mean matrices. For each mean matrix, it generates `cn^2` Gaussian random matrices (with identity covariance matrix), and compute the empirical mean and empirical covariance as consistent estimators. It then combines those estimates into a list of sub-lists that include the estimated mean `mu.bar` and covariance `cov.bar`.

7. For more details, please refer to the vignette: `browseVignettes('eigTest')` and the paper 'Testing Simultaneous Diagonalizability'.

## Code
The folder [`code`](code) includes the scripts for different simulation studies and applications in the paper.

In particular, the simulation studies include examples for commutator-based test, LLR test, multi-sample eigenvector test and partial test. The exploratory script for two-sample tests under relatively high-dimensional setting is also provided. Applications include VAR examples on macroeconomic data of 8 countries, and Markovian model on Hudson river discharge data.

## Data
The folder [`data`](data) includes the `R` objects of:
1. the array of estimated VAR coefficient matrices for 8 countries `countryCoeff.RData`;
2. the estimated covariance of VAR coefficient matrices for 8 countries `countryCovar.RData`;
3. the weekly streamflow discharge data of Hudson river `hudsonWeekly.RData`;
4. the daily streamflow discharge data of Hudson river `hudsonDaily.RData`.
5. the replicates of raw p-values from simulation studies for test power analysis and histogram visualizations.

The first 4 datasets are also available from the package `eigTest`.

The only `csv` file [`M2_KOR_monthly_SA.csv`](data/M2_KOR_monthly_SA.csv) in the folder has the Korean money supply data, which cannot be directly read by functions in `R` from the Bank of Korea website.

## Output
The folder [`output`](output) includes the replicates of raw p-values from simulation studies in folder [`code`](code), and the intermediate output of macroeconomic multivariate time series data collected and processed by [`VAR_application.R`](code/VAR_application.R). The sub-folder [`Plots`](output/Plots) includes generated plots from both simulations and applications in `code`, with either [`tikz`](output/Plots/tikz) or [`png`](output/Plots/png) format.
