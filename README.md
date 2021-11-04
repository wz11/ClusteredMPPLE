
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ClusteredMPPLE

<!-- badges: start -->

<!-- badges: end -->

The ClusteredMPPLE conducts semiparametric marginal regression for
clustered competing risks data with missing cause of failure. It
implements the proposed maximum partial pseudolikelihood estimation
method by Zhou, W. et al.(2021) based on semiparametric marginal
proportional cause-specific hazards model. The standard errors for the
estimated regression coefficients are calculated based on the close-form
variance or bootstrapping. It also calculates and generates cumulative
incidence functions and cumulative residual processes plots.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wz11/ClusteredMPPLE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(survival)
library(ClusteredMPPLE)

# Load simulated dataset
data("simulated_data")

# Fit the model 
fit <- ccr_smreg(data=simulated_data, formula11 =  y ~ x + Z1 + Z2, formula12 =  y ~ x + Z1 + Z2, formula21 = Surv(x, d) ~ Z1 + Z2, formula22 = Surv(x, d) ~ Z1 + Z2, w = TRUE, var.method = "BS", nboot = 1000)

#> please wait...

# Obtain model summary
summary(fit)
#> 
#> cause 1: 
#>          beta         se      p-value exp(beta) 95% CI lower limit
#> Z1 -0.1949448 0.02551021 2.141514e-14 0.8228801          0.7827478
#> Z2 -0.0237731 0.06768892 7.254300e-01 0.9765073          0.8551799
#>    95% CI upper limit
#> Z1          0.8650701
#> Z2          1.1150477
#> 
#> cause 2: 
#>           beta         se     p-value exp(beta) 95% CI lower limit
#> Z1 -0.05631217 0.03780163 0.136309589 0.9452440          0.8777413
#> Z2 -0.36750156 0.13069940 0.004926395 0.6924622          0.5359725
#>    95% CI upper limit
#> Z1          1.0179380
#> Z2          0.8946428

# Calculate and plot cumulative incidence function
pfit <- cif.ccr_smreg(fit, band = TRUE, sims = 1000)
plot(pfit)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

# Calculate and plot cumulative residual process
rfit <- crp.ccr_smreg(fit)
plot(rfit)
```

<img src="man/figures/README-example-2.png" width="100%" />

## Reference

  - Zhou, W., Bakoyannis, G., Zhang, Y., & Yiannoutsos, C. T. (2021).
    Semiparametric Marginal Regression for Clustered Competing Risks
    Data with Missing Cause of Failure. arXiv preprint arXiv:2104.09090.
