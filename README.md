# SDAtools: A toolkit for analysing SDA decompositions

This is a user friendly R package for facilitating the use of [SDA](https://jmarchini.org/sda/), including data formatting, importing results and analysis.

## Installation

```R
# install.packages("remotes")
remotes::install_github("marchinilab/SDAtools")
```

## Quick Example Usage
```R
library(SDAtools)

run_SDA(out = "simulation_results",
        data = "simulated.data",
        N = 100)

results <- load_results(results_folder = "simulation_results", iteration = 5000)

check_convergence(results)

highest_genes(results, component = 1)
```

For a full guide on how to use this package please see the [vignette](vignettes/vignette.md).

## Directory Organisation
The [R/](R/) directory contains source code of the functions.

The [man/](man/) directory contains the manual pages for the functions, compiled by roxygen.

The [vignettes/](vignettes/) directory contains tutorial exemplifying usage.

The [tests/](tests/) directory contains the unit tests which are carried out by the testthat R package.

Contributions of any size or form are welcome!