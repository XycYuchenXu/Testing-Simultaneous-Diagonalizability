This `output` folder includes output objects from manuscript implementations with multiple formats. In particular:
1. `RData`:

    a. `commutatorTest.RData`: raw p-values for commutator-based test simulations (Section 3 & 6.2);
    
    b. `LLRTest.RData`: raw p-values for log-likelihood ratio test simulations (Section B at supplementary);
    
    c. `highD_Test.RData`: raw p-values for two-sample test simulations under high-dimensional setting (with increasing and relatively larger dimensions, Section A.3 at supplementary);
    
    d. `multiTest.RData`: raw p-values for multi-sample test simulations (Section 4 & 6.3);
    
    e. `partialTest.RData`: raw p-values for partial test simulations (Section 5 & 6.4);
    
    f. `countryMacro.RData`: the `list` of macroeconomic quarterly time series for 8 countries - the intermediate output after online retrieval and initial processing.
    
    These first 5 files of p-values are also viewed as data to generate summary histograms and error tables, and hence having duplicates in the `data` folder.
2. `Error Tables.xlsx`: the summarised type I/II errors for all simulations setups;
3. `Plots`:

    a. `tikz`: folder of source `tikz` files for compiling the manuscript;
    
    b. `png`: folder of `png` files for easy visualization.
