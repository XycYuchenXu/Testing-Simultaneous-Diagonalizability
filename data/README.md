This `data` folder includes:

1. `commutatorTest.RData`: raw p-values for commutator-based test simulations (Section 3 & 6.2);

2. `LLRTest.RData`: raw p-values for log-likelihood ratio test simulations (Section A at supplementary);

3. `highD_Test.RData`: raw p-values for two-sample test simulations under high-dimensional setting (with increasing and relatively larger dimensions, Section B.3 at supplementary);

4. `multiTest.RData`: raw p-values for multi-sample test simulations (Section 4 & 6.3);

5. `partialTest.RData`: raw p-values for partial test simulations (Section 5 & 6.4);

6. `M2_KOR_monthly_SA.csv`: the publicly available seasonally adjusted monthly money aggregate M2 for Korea, accessed from Bank of Korea.

7. `countryCoeff.RData`: the array of estimated VAR coefficient matrices for 8 countries;

8. `countryCovar.RData`: the estimated covariance of VAR coefficient matrices for 8 countries;

9. `hudsonWeekly.RData`: the weekly streamflow discharge data of Hudson river;

10. `hudsonDaily.RData`: the daily streamflow discharge data of Hudson river;

The first 5 datasets may also be considered as output of the project and hence have duplicates in the [output](../output) folder. The M2 data of Korea may not be directly URL-accessed by `R` from the website hence this separate `csv` table is included. The latter 4 datasets are also available from the `R` package [`eigTest`](https://github.com/XycYuchenXu/eigTest).
