
<!-- README.md is generated from README.Rmd. Please edit this file -->

# gddlm – Robust double linear models using gamma divergence

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

Please check the help file for the main `gddlm` function after
installation. Alternatively, you can follow the basic example below.

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
set.seed(012026)
n <- 250
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

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

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
```

``` r

fit_gddglm_choosegamma <- gddlm(mu_formula = make_mu_formula,
                                sigma_formula = make_sigma_formula,
                                data = data.frame(response = Y, X[,-1]),
                                gamma_tuning = runcv$gamma_optimal)
```

Demonstration of some of the functions available in the package

``` r
fit_gddglm_choosegamma %>% 
    str
#> List of 15
#>  $ call                : language gddlm(mu_formula = make_mu_formula, sigma_formula = make_sigma_formula,      data = data.frame(response = Y, X[, | __truncated__
#>  $ mu_formula          :Class 'formula'  language response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#>   .. ..- attr(*, ".Environment")=<environment: 0x5c80cbde12e8> 
#>  $ sigma_formula       :Class 'formula'  language ~X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#>   .. ..- attr(*, ".Environment")=<environment: 0x5c80cbe95960> 
#>  $ model               :'data.frame':    250 obs. of  11 variables:
#>   ..$ response: num [1:250] 0.02842 -1.13633 0.26403 1.51584 -0.00728 ...
#>   ..$ X1      : num [1:250] -0.743 -1.604 -0.867 -0.973 -0.34 ...
#>   ..$ X2      : num [1:250] -0.383 -1.58 0.505 -0.894 -0.606 ...
#>   ..$ X3      : num [1:250] 0.08 -0.3233 0.2707 -0.6501 0.0689 ...
#>   ..$ X4      : num [1:250] 0.54 -0.198 0.597 -0.684 -1.314 ...
#>   ..$ X5      : num [1:250] 0.782 -1.444 0.765 -0.977 -0.963 ...
#>   ..$ X6      : num [1:250] 0.0218 -1.0271 0.6141 -0.6685 0.1256 ...
#>   ..$ X7      : num [1:250] -1.2093 -0.3664 0.3461 -0.132 0.0804 ...
#>   ..$ X8      : num [1:250] -0.987 -0.552 -0.296 -0.351 0.268 ...
#>   ..$ X9      : num [1:250] -0.0254 -0.8484 0.4772 -0.6726 0.5565 ...
#>   ..$ X10     : num [1:250] 0.294 -1.241 0.768 -0.739 1.281 ...
#>   ..- attr(*, "terms")=Classes 'terms', 'formula'  language response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#>   .. .. ..- attr(*, "variables")= language list(response, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
#>   .. .. ..- attr(*, "factors")= int [1:11, 1:10] 0 1 0 0 0 0 0 0 0 0 ...
#>   .. .. .. ..- attr(*, "dimnames")=List of 2
#>   .. .. .. .. ..$ : chr [1:11] "response" "X1" "X2" "X3" ...
#>   .. .. .. .. ..$ : chr [1:10] "X1" "X2" "X3" "X4" ...
#>   .. .. ..- attr(*, "term.labels")= chr [1:10] "X1" "X2" "X3" "X4" ...
#>   .. .. ..- attr(*, "order")= int [1:10] 1 1 1 1 1 1 1 1 1 1
#>   .. .. ..- attr(*, "intercept")= int 1
#>   .. .. ..- attr(*, "response")= int 1
#>   .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>   .. .. ..- attr(*, "predvars")= language list(response, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10)
#>   .. .. ..- attr(*, "dataClasses")= Named chr [1:11] "numeric" "numeric" "numeric" "numeric" ...
#>   .. .. .. ..- attr(*, "names")= chr [1:11] "response" "X1" "X2" "X3" ...
#>  $ mu_coefficients     : Named num [1:11] -1.045 -1.11 1.253 -0.656 0.441 ...
#>   ..- attr(*, "names")= chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ sigma_coefficients  : Named num [1:11] 0.474 0.484 -0.144 0.156 -0.478 ...
#>   ..- attr(*, "names")= chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ mu_fitted           : num [1:250] -0.149 -1.551 0.98 -1.276 -2.366 ...
#>  $ sigma_fitted        : num [1:250] 0.863 1.334 0.757 1.69 1.914 ...
#>  $ raw_residuals       : num [1:250] 0.177 0.415 -0.716 2.792 2.359 ...
#>  $ pearson_residuals   : num [1:250] 0.205 0.311 -0.946 1.652 1.232 ...
#>  $ final_weights       : num [1:250] 1.389 1.209 1.272 0.759 0.876 ...
#>  $ gamma_tuning        : num 0.3
#>  $ gamma_divergence    : num -446
#>  $ V_mu_coefficients   : num [1:11, 1:11] 0.01909 0.01114 -0.00399 0.00179 -0.00784 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  $ V_sigma_coefficients: num [1:11, 1:11] 0.003675 0.000201 0.000582 -0.000734 0.000143 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>   .. ..$ : chr [1:11] "(Intercept)" "X1" "X2" "X3" ...
#>  - attr(*, "class")= chr "gddlm"
```

``` r

summary(fit_gddglm_choosegamma)
#> Double linear model fitted using gamma divergence with gamma = 0.3 
#> Formula:
#> response ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#> <environment: 0x5c80cbde12e8>
#> ~X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
#> <environment: 0x5c80cbe95960>
#> 
#> Parametric coefficients (mu and sigma):
#>                    Estimate Std. Error z value Pr(>|z|)    
#> mu_(Intercept)    -1.045292   0.138169  -7.565 3.87e-14 ***
#> mu_X1             -1.110320   0.145234  -7.645 2.09e-14 ***
#> mu_X2              1.252945   0.190682   6.571 5.00e-11 ***
#> mu_X3             -0.655773   0.174682  -3.754 0.000174 ***
#> mu_X4              0.441316   0.156681   2.817 0.004853 ** 
#> mu_X5              0.335970   0.185683   1.809 0.070393 .  
#> mu_X6             -0.097287   0.172607  -0.564 0.573007    
#> mu_X7              0.056360   0.180755   0.312 0.755192    
#> mu_X8             -0.190879   0.155283  -1.229 0.218984    
#> mu_X9              0.208394   0.139987   1.489 0.136574    
#> mu_X10            -0.037177   0.130372  -0.285 0.775522    
#> sigma_(Intercept)  0.473962   0.060622   7.818 5.35e-15 ***
#> sigma_X1           0.484237   0.091212   5.309 1.10e-07 ***
#> sigma_X2          -0.144390   0.120539  -1.198 0.230966    
#> sigma_X3           0.155572   0.119813   1.298 0.194130    
#> sigma_X4          -0.478307   0.112075  -4.268 1.97e-05 ***
#> sigma_X5           0.008951   0.137674   0.065 0.948159    
#> sigma_X6          -0.026517   0.127788  -0.208 0.835611    
#> sigma_X7           0.187717   0.116311   1.614 0.106546    
#> sigma_X8          -0.267821   0.124803  -2.146 0.031878 *  
#> sigma_X9           0.269928   0.129544   2.084 0.037189 *  
#> sigma_X10         -0.365758   0.116662  -3.135 0.001717 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> 
#> Residual degrees of freedom: 228 
#> Gamma divergence: -446.0563
```

``` r

appraise(fit_gddglm_choosegamma)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

``` r

coef(fit_gddglm_choosegamma)
#> (Intercept)          X1          X2          X3          X4          X5 
#> -1.04529195 -1.11031958  1.25294484 -0.65577341  0.44131644  0.33597042 
#>          X6          X7          X8          X9         X10 
#> -0.09728660  0.05635957 -0.19087916  0.20839436 -0.03717663
```

``` r
coef(fit_gddglm_choosegamma, type = "sigma")
#>  (Intercept)           X1           X2           X3           X4           X5 
#>  0.473962012  0.484237441 -0.144390181  0.155571601 -0.478306825  0.008951375 
#>           X6           X7           X8           X9          X10 
#> -0.026517334  0.187716786 -0.267820596  0.269927990 -0.365758004
```

``` r

simulate(fit_gddglm_choosegamma, nsim = 10) %>% 
    str
#> 'data.frame':    250 obs. of  10 variables:
#>  $ sim_1 : num  -1.359 -0.207 1.076 -2.468 -2.218 ...
#>  $ sim_2 : num  -0.627 -1.939 1.142 -0.408 -3.777 ...
#>  $ sim_3 : num  -0.229 0.449 0.133 -1.167 -2.497 ...
#>  $ sim_4 : num  1.278 0.164 0.852 -1.744 -2.32 ...
#>  $ sim_5 : num  -0.26 -1.417 0.713 -0.631 -1.829 ...
#>  $ sim_6 : num  0.444 -3.855 1.372 2.247 -3.665 ...
#>  $ sim_7 : num  0.294 -1.776 2.388 -1.877 -5.421 ...
#>  $ sim_8 : num  0.0443 -2.8441 -0.0128 -3.489 -5.0201 ...
#>  $ sim_9 : num  -0.0269 -4.4495 1.3718 -0.8503 -0.8624 ...
#>  $ sim_10: num  0.382 1.092 1.636 -1.481 -1.743 ...
#>  - attr(*, "seed")= int [1:626] 10403 40 -1053275711 1919596251 509790488 1295980936 -960296891 94531324 -380244184 -532208171 ...
```

``` r

predict(fit_gddglm_choosegamma, predictor = "mu", se.fit = TRUE) %>% 
    head
#>          fit    se.fit      lower       upper
#> 1 -0.1489594 0.2277376 -0.5953169  0.29739812
#> 2 -1.5512057 0.3037957 -2.1466344 -0.95577696
#> 3  0.9795348 0.2255716  0.5374227  1.42164695
#> 4 -1.2763706 0.1910713 -1.6508635 -0.90187764
#> 5 -2.3659976 0.2775728 -2.9100302 -1.82196486
#> 6 -0.7416024 0.4070302 -1.5393669  0.05616203
```

``` r
predict(fit_gddglm_choosegamma, predictor = "sigma", se.fit = TRUE) %>% 
    head
#>         fit   se.fit     lower    upper
#> 1 0.8632558 1.200745 0.6031408 1.235550
#> 2 1.3339486 1.215766 0.9095707 1.956328
#> 3 0.7565296 1.202796 0.5268083 1.086424
#> 4 1.6898861 1.116319 1.3620536 2.096625
#> 5 1.9137916 1.235816 1.2637721 2.898148
#> 6 3.6122831 1.287683 2.2006950 5.929304
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
