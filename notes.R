rm(list = ls())
library(tidyverse)
library(GJRM)
library(mgcv)
library(mvtnorm)

# Simulation setting and data generation
set.seed(012026)
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
om <- 0.1
out <- rbinom(n, 1, om)
Y[out==1] <- Y[out==1] + 50

hist(Y, breaks = 20)


# Fit gamma divergence double linear models
make_mu_formula <- paste0("response ~" , paste0("X", 1:p, collapse = " + ")) %>% formula
make_sigma_formula <- paste0("~" , paste0("X", 1:p, collapse = " + ")) %>% formula


# fit_nonrobust <- gam(list(make_mu_formula, make_sigma_formula),
#                      data = data.frame(response = Y, X[,-1]),
#                      family = gaulss())
#
# fit_gddglm <- gddlm(mu_formula = make_mu_formula,
#                     sigma_formula = make_sigma_formula,
#                     data = data.frame(response = Y, X[,-1]),
#                     gamma_tuning = 0.1)
#
#
# fit_nonrobust$coefficients
# summary(fit_gddglm)


runcv <- cv_gddlm(mu_formula = make_mu_formula,
                  sigma_formula = make_sigma_formula,
                  data = data.frame(response = Y, X[,-1]))

fit_gddglm_choosegamma <- gddlm(mu_formula = make_mu_formula,
                                sigma_formula = make_sigma_formula,
                                data = data.frame(response = Y, X[,-1]),
                                gamma_tuning = runcv$gamma_optimal)

fit_gddlm_lsa <- gddlm_lsa(object = fit_gddglm_choosegamma)
apply(fit_gddlm_lsa$information, 2, which.min)
fit_gddlm_lsa$mu_coefficients_path[,which.min(fit_gddlm_lsa$information$BIC)]
fit_gddlm_lsa$sigma_coefficients_path[,which.min(fit_gddlm_lsa$information$BIC)]


fit_gddglm_gamma0 <- gddlm(mu_formula = make_mu_formula,
                                sigma_formula = make_sigma_formula,
                                data = data.frame(response = Y, X[,-1]),
                                gamma_tuning = 1e-12)

fit_gddlm_gamma0_lsa <- gddlm_lsa(object = fit_gddglm_gamma0)
apply(fit_gddlm_gamma0_lsa$information, 2, which.min)
fit_gddlm_gamma0_lsa$mu_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)]
fit_gddlm_gamma0_lsa$sigma_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)]




function() {
    mu_formula = make_mu_formula
    sigma_formula = make_sigma_formula
    data = data.frame(response = Y, X[,-1])
    gamma_tuning = 1e-12
    mu_bar_lambda = 0.2
    sigma_bar_lambda = 0.2
    control = list(max_iteration = 1e3, tol = 1e-5, dampening = c(1,1), trace = FALSE)


    object = fit_gddglm_choosegamma
    nlambda = 100
    lambda_min_ratio = 1e-6
    lambda = NULL
    min_df = 0
    num_cores = NULL
    prior_weights = NULL
    control = list(maxit = 1000, eps = 1e-6, round_eps = 1e-6)
    }
