# Testing Simultaneous Diagonalizability

This is the reproducible repository for the paper. It mainly relys on the `R` package `eigTest` developed by the authors.

## Setup
The `R` package `eigTest` is available at Github for testing simultaneous diagonalizability.

To install, run the following code in `R`:
`devtools::install_github('XycYuchenXu/eigTest', build_vignettes = T)`

## Usage
That package has the following functionalities:

1. Check whether the means of two random square matrices (`matrixA` and `matrixB`) with dimension `d`-by-`d` are simultaneously diagonalizable, either considering the commutator of the two matrices (use function `commutatorTest`), or using the log-likelihood ratio test framework (`projTest`), given asymptotic limiting covariance matrices (`covMatA` and `covMatB`) and a convergence rate (`cn`). For both functions, one needs to input an array of two matrices (`A` such that `A[1,,] = matrixA` and `A[2,,] = matrixB`), an array of two limiting covariance matrices (`cov.arr` such that `cov.arr[1,,] = covMatA` and `cov.arr[2,,] = covMatB`), and the convergence rate (`cn`).

2. Given an array of square random matrices (`A`), estimate a common eigenvector matrix `V` (by supplying the array of matrices `A` for function `JDTE`).

3. Test whether the array of means matrices (`A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) for random matrices are simultaneously diagonalizable, by checking whether the eigenvector matrix estimated above (`V`) is indeed the common eigenvector matrix. Use function `eigTest` with input of two arrays (`A` and `cov.arr`) and the convergence rate `cn`.

4. Given an array of square random matrices (`A` with dimension `p`-by-`d`-by-`d`), estimate `k` (`k` ≤ `d`) common eigenvectors `Vk` (by supplying the array of matrices `A` and parameter `k` for function `expmPartSchur`).

5. Test whether the array of the random estimates (`A` with dimension `p`-by-`d`-by-`d` and `p` ≥ 2) for mean matrices share `k` of the all `d` eigenvectors (`k` ≤ `d`), by checking whether the `k` estimated common eigenvectors `Vk` from optimization above are indeed the common eigenvectors. Use function `partialTest` with input of two arrays (`A` and `cov.arr`), the convergence rate `cn`, and number of common eigenvectors `k`. `Vk` is optional here as `partialTest` can call the optimization function.

6. Generate Gaussian samples for simulations. The function `generateMeans` is used to generate `p` mean matrices with dimension `d`-by-`d` and noise levels `snr` of common eigenvector `V`. Function `simuSamples` takes the output from `generateMeans` as input for mean matrices. For each mean matrix, it generates `cn^2` Gaussian random matrices (with identity covariance matrix), and compute the sample mean and the sample covariance matrix as consistent estimators. It then combines those estimates into a list of sub-lists that include estimated mean `mu.bar` and covariance `cov.bar`.

7. For more details, read the vignette: `browseVignettes('eigTest')` and the paper 'Testing Simultaneous Diagonalizability'.

## Code
The folder `code` includes the scripts for different simulation studies and applications in the paper.

In particular, the simulation studies include examples for commutator-based test, LLR test, multi-sample eigenvector test and partial test. And applications include VAR examples on macroeconomic data of 8 countries and Markovian model on Hudson river discharge data.

## Data
The folder `data` includes the `R` objects of:
1. the macroeconomic time series of 8 countries `countryMacro.RData`;
2. the array of estimated VAR coefficient matrices for 8 countries `countryCoeff.RData`;
3. the estimated covariance of VAR coefficient matrices for 8 countries `countryCovar.RData`;
4. the weekly streamflow discharge data of Hudson river `hudsonWeekly.RData`;
5. the daily streamflow discharge data of Hudson river `hudsonDaily.RData`.
These datasets are also available from the package `eigTest`.

## Output
The folder `output` includes the replicates of raw p-values from simulation studies in folder `code`. The sub-folder `Plots` includes generated plots from both simulations and applications in `code`, with either `tikz` or `png` format.
