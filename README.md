# vetperf

<!-- badges: start -->
<!-- badges: end -->

The goal of vetperf is to provide functions for extracting perfusion parameters from
measurements obtained from veterinary MRI images.

## Installation

You can install the released version of vetperf from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("vetperf")
```

## Example

Load your data from csv (check out the great data.table package!) and run vetperf.data on
your data. You will get a table holding the perfusion parameters.

``` r
library(data.table)
library(vetperf)
animals <- fread('path/to/your/data.csv', stringsAsFactors = TRUE)
perfdata = vetperf.data(animals, arrival_frame = 9, delta_t = 1.6, echo_time = 0.03, baseroi = "cs")
```
## More information

Can be found in the package vignette. For an explanation on the format of the data frame, take a 
look at the examples or at the code documentation.
