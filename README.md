
<!-- README.md is generated from README.Rmd. Please edit this file -->

# gddlm: Robust double linear models using gamma divergence

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)

<!-- badges: end -->

`gddlm` is an `R` package associated with the manuscript “Robust double
linear models using gamma divergence” by
[Hui](https://francishui.netlify.app/), and
[Sugusawa](https://sites.google.com/view/ssugasawa/), which is currently
in preparation.

# Installation

Currently, `gddlm` is available and can be installed from Github with
the help of `pak` package using:

    pak::pkg_install("fhui28/gddlm")

<!-- Alternatively, or if the above does not work, you may download a (hopefully!) stable release of `COQUE` by choosing the latest release on the right hand side of this Github webpage. -->

# Getting started and example usage

Please check the help files for the main `gddlm` and `gddlm_lsa`
functions after installation. Alternatively, you can follow the basic
example below.

``` r
library(tidyverse)
library(mgcv)
library(mvtnorm)
library(gratia)
library(gddlm)
```

Simulate data from a double linear model with (average 5%) outliers in
the response variable

``` r
set.seed(112026)
n <- 500
p <- 10
H <- abs(outer(1:p, 1:p, "-"))
H <- 0.75^H
X <- cbind(1, rmvnorm(n, sigma = H))
Beta_true <- c(-1, -1, 1, -0.5, 0.5, 0.2, -0.2, rep(0, p - 6))
Delta_true <- c(0.5, 0.5, 0, 0, -0.5, 0, 0, 0.25, -0.25, rep(0, p - 8))

## Generate responses
mu_true <- as.vector(X %*% Beta_true)
sig_true <- exp(as.vector(X %*% Delta_true))
Y <- rnorm(n, mean = mu_true, sd = sig_true)

## Create outliers
om <- 0.05
out <- rbinom(n, 1, om)
Y[out==1] <- Y[out==1] + 50

hist(Y, breaks = 20)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" alt="" width="100%" />

Fit gamma divergence double linear models, select the gamma tuning
parameter via cross-validation

``` r
make_mu_formula <- paste0("response ~" , paste0("X", 1:p, collapse = " + ")) %>%
    formula
make_sigma_formula <- paste0("~" , paste0("X", 1:p, collapse = " + ")) %>%
    formula

runcv <- cv_gddlm(mu_formula = make_mu_formula,
                  sigma_formula = make_sigma_formula,
                  data = data.frame(response = Y, X[,-1]))
#> Onto fold 1...
#> Onto fold 2...
#> Onto fold 3...
#> Onto fold 4...
#> Onto fold 5...
#> Onto fold 6...
#> Onto fold 7...
#> Onto fold 8...
#> Onto fold 9...
#> Onto fold 10...

fit_gddglm_choosegamma <- gddlm(mu_formula = make_mu_formula,
                                sigma_formula = make_sigma_formula,
                                data = data.frame(response = Y, X[,-1]),
                                gamma_tuning = runcv$gamma_optimal_mean)
```

Demonstration of some of the functions available in the package

``` r
fit_gddglm_choosegamma %>% 
    str
#> List of 16
#>  $ call                : language gddlm(mu_formula = make_mu_formula, sigma_formula = make_sigma_formula,      data = data.frame(response = Y, X[, | __truncated__
#>  $ mu_formula          :Class 'formula'  language response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#>   .. ..- attr(*, ".Environment")=<environment: 0x5f8ebbddca80> 
#>  $ sigma_formula       :Class 'formula'  language ~X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#>   .. ..- attr(*, ".Environment")=<environment: 0x5f8ebbe0c708> 
#>  $ model               :'data.frame':    500 obs. of  11 variables:
#>   ..$ response: num [1:500] -2.361 -2.709 -0.932 -0.971 -1.028 ...
#>   ..$ X1      : num [1:500] 1.2748 1.0061 0.0843 0.3459 -0.8023 ...
#>   ..$ X2      : num [1:500] 1.007 0.773 -0.757 -0.763 0.522 ...
#>   ..$ X3      : num [1:500] 0.92 1.483 -0.559 -1.259 1.073 ...
#>   ..$ X4      : num [1:500] 1.596 1.324 -1.009 -0.748 -0.457 ...
#>   ..$ X5      : num [1:500] 1.0546 1.52736 0.00438 -0.37277 -0.53995 ...
#>   ..$ X6      : num [1:500] 0.178 1.217 -0.401 -1.178 -0.904 ...
#>   ..$ X7      : num [1:500] 0.285 0.66 0.108 -1.695 -0.124 ...
#>   ..$ X8      : num [1:500] -0.283 0.155 0.581 -2.436 0.285 ...
#>   ..$ X9      : num [1:500] 0.953 -0.427 0.431 -1.433 0.776 ...
#>   ..$ X10     : num [1:500] 0.802 1.126 1.378 -1.335 0.806 ...
#>  $ mu_coefficients     : Named num [1:11] -1.052 -1.04 1.134 -0.83 0.654 ...
#>   ..- attr(*, "names")= chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ sigma_coefficients  : Named num [1:11] 0.421 0.436 0.228 -0.175 -0.567 ...
#>   ..- attr(*, "names")= chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ mu_fitted           : num [1:500] -0.857 -1.459 -2.023 -1.614 -0.727 ...
#>  $ sigma_fitted        : num [1:500] 1.38 1.18 2.24 3.27 1.14 ...
#>  $ raw_residuals       : num [1:500] -1.504 -1.25 1.091 0.644 -0.301 ...
#>  $ pearson_residuals   : num [1:500] -1.089 -1.056 0.487 0.197 -0.263 ...
#>  $ final_weights       : num [1:500] 1.003 1.082 0.999 0.894 1.351 ...
#>  $ gamma_tuning        : num 0.4
#>  $ gamma_divergence    : num -871
#>  $ V_full              : num [1:22, 1:22] 0.009135 0.00229 0.001004 -0.000225 -0.006487 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:22] "(Intercept)" "X1" "X2" "X3" ...
#>   .. ..$ : chr [1:22] "(Intercept)" "X1" "X2" "X3" ...
#>  $ V_mu_coefficients   : num [1:11, 1:11] 0.009135 0.00229 0.001004 -0.000225 -0.006487 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ V_sigma_coefficients: num [1:11, 1:11] 1.80e-03 1.02e-03 -2.04e-03 1.96e-03 4.92e-05 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  - attr(*, "class")= chr "gddlm"

summary(fit_gddglm_choosegamma)
#> Double linear model fitted using gamma divergence with gamma = 0.4 
#> Formula:
#> response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#> <environment: 0x5f8ebbddca80>
#> ~X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#> <environment: 0x5f8ebbe0c708>
#> 
#> Parametric coefficients (mu and sigma):
#>                    Estimate Std. Error z value Pr(>|z|)    
#> mu_(Intercept)    -1.051772   0.095577 -11.004  < 2e-16 ***
#> mu_X1             -1.039828   0.087575 -11.874  < 2e-16 ***
#> mu_X2              1.133576   0.118410   9.573  < 2e-16 ***
#> mu_X3             -0.830221   0.124197  -6.685 2.31e-11 ***
#> mu_X4              0.653921   0.142897   4.576 4.74e-06 ***
#> mu_X5              0.243204   0.106373   2.286 0.022235 *  
#> mu_X6             -0.243222   0.169397  -1.436 0.151055    
#> mu_X7             -0.069251   0.112345  -0.616 0.537620    
#> mu_X8              0.120493   0.121449   0.992 0.321136    
#> mu_X9             -0.093130   0.104630  -0.890 0.373418    
#> mu_X10             0.035339   0.082598   0.428 0.668769    
#> sigma_(Intercept)  0.420895   0.042391   9.929  < 2e-16 ***
#> sigma_X1           0.436283   0.080195   5.440 5.32e-08 ***
#> sigma_X2           0.228361   0.114474   1.995 0.046056 *  
#> sigma_X3          -0.174513   0.111419  -1.566 0.117282    
#> sigma_X4          -0.567367   0.086588  -6.552 5.66e-11 ***
#> sigma_X5           0.074916   0.087645   0.855 0.392683    
#> sigma_X6          -0.067892   0.106915  -0.635 0.525422    
#> sigma_X7           0.320068   0.088858   3.602 0.000316 ***
#> sigma_X8          -0.231268   0.084389  -2.740 0.006135 ** 
#> sigma_X9           0.003943   0.085225   0.046 0.963101    
#> sigma_X10         -0.056193   0.070130  -0.801 0.422976    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> Residual degrees of freedom: 478 
#> Gamma divergence: -870.7739

appraise(fit_gddglm_choosegamma)
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the gddlm package.
#>   Please report the issue at <https://github.com/fhui28/gddlm/issues>.
#> This warning is displayed once per session.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
```

<img src="man/figures/README-unnamed-chunk-5-1.png" alt="" width="100%" />

``` r

coef(fit_gddglm_choosegamma)
#> (Intercept)          X1          X2          X3          X4          X5 
#> -1.05177231 -1.03982829  1.13357610 -0.83022127  0.65392125  0.24320395 
#>          X6          X7          X8          X9         X10 
#> -0.24322213 -0.06925101  0.12049300 -0.09313003  0.03533871
coef(fit_gddglm_choosegamma, type = "sigma")
#>  (Intercept)           X1           X2           X3           X4           X5 
#>  0.420894849  0.436283286  0.228360716 -0.174513315 -0.567367015  0.074915934 
#>           X6           X7           X8           X9          X10 
#> -0.067892115  0.320068324 -0.231268126  0.003942675 -0.056193169

simulate(fit_gddglm_choosegamma, nsim = 10) %>% 
    str
#> 'data.frame':    500 obs. of  10 variables:
#>  $ sim_1 : num  -1.75 -1.41 -5.39 -3.94 -1.06 ...
#>  $ sim_2 : num  -0.00866 -1.88531 0.05362 -1.70807 0.70808 ...
#>  $ sim_3 : num  -1.324 -0.814 0.889 -8.047 0.414 ...
#>  $ sim_4 : num  -3.3469 0.0211 0.7229 -0.9174 -2.66 ...
#>  $ sim_5 : num  -0.2126 -1.7837 -0.0696 5.2576 -1.5093 ...
#>  $ sim_6 : num  -1.676 -1.164 -2.434 -5.23 0.268 ...
#>  $ sim_7 : num  -1.847 -0.466 -3.129 0.267 -1.458 ...
#>  $ sim_8 : num  -1.598 -1.397 -6.67 -0.339 -1.259 ...
#>  $ sim_9 : num  -1.68 -3.215 0.103 -0.344 0.544 ...
#>  $ sim_10: num  -0.76 -0.559 -1.284 -3.502 0.681 ...
#>  - attr(*, "seed")= int [1:626] 10403 80 638472023 -924926401 -1111118522 1792087268 1698227145 280349103 1720174138 -959805469 ...

predict(fit_gddglm_choosegamma, predictor = "mu", se.fit = TRUE) %>% 
    head
#>          fit    se.fit     lower      upper
#> 1 -0.8565576 0.2455725 -1.337871 -0.3752443
#> 2 -1.4589825 0.1755096 -1.802975 -1.1149899
#> 3 -2.0233802 0.1683308 -2.353303 -1.6934579
#> 4 -1.6141223 0.2525456 -2.109103 -1.1191420
#> 5 -0.7272791 0.1996628 -1.118611 -0.3359471
#> 6 -0.2980597 0.2295237 -0.747918  0.1517985
predict(fit_gddglm_choosegamma, predictor = "sigma", se.fit = TRUE) %>% 
    head
#>        fit   se.fit     lower    upper
#> 1 1.381538 1.178141 1.0018875 1.905052
#> 2 1.183412 1.179174 0.8567342 1.634654
#> 3 2.241801 1.136381 1.7449065 2.880197
#> 4 3.267212 1.148173 2.4920964 4.283412
#> 5 1.144709 1.200304 0.8003628 1.637205
#> 6 1.534821 1.175897 1.1172150 2.108526
```

# If you find any bugs and issues…

If you find something that looks like a bug/issue, please use Github
issues and post it up there. As much as possible, please include in the
issue:

1.  A description of the bug/issue;
2.  Paste-able code along with some comments that reproduces the problem
    e.g., using the
    [reprex](https://cran.r-project.org/web/packages/reprex/index.html)
    package. If you also have an idea of how to fix the problem, then
    that is also much appreciated.
3.  Required data files etc…

Alternatively, please contact the corresponding author at
<francis.hui@anu.edu.au>
