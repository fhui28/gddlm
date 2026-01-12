rm(list = ls())
library(tidyverse)
library(GJRM)
library(mgcv)
library(mvtnorm)
library(gamlss.lasso)
library(DoubleRobGam) #remotes::install_github("ilapros/DoubleRobGam")


# Simulation setting and data generation
set.seed(012026)
n <- 500
p <- 5*ceiling(n^0.3)
H <- abs(outer(1:p, 1:p, "-"))
H <- 0.5^H
X <- cbind(1, rmvnorm(n, sigma = H))
Beta_true <- c(-1, -1, 1, -0.5, 0.5, 0, 0, 0.25, rep(0, p - 7))
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


# Fit gamma divergence double linear models along with several other methods
make_mu_formula <- paste0("response ~" , paste0("X", 1:p, collapse = " + ")) %>% formula
make_sigma_formula <- paste0("~" , paste0("X", 1:p, collapse = " + ")) %>% formula

fit_lmls <- gam(list(make_mu_formula, make_sigma_formula),
                data = data.frame(response = Y, X[,-1]),
                family = gaulss(),
                method = "ML")


runcv <- cv_gddlm(mu_formula = make_mu_formula,
                  sigma_formula = make_sigma_formula,
                  data = data.frame(response = Y, X[,-1]))

fit_gddlm_choosegamma <- gddlm(mu_formula = make_mu_formula,
                                sigma_formula = make_sigma_formula,
                                data = data.frame(response = Y, X[,-1]),
                                gamma_tuning = runcv$gamma_optimal)


fit_gddlm0 <- gddlm(mu_formula = make_mu_formula,
                    sigma_formula = make_sigma_formula,
                    data = data.frame(response = Y, X[,-1]),
                    gamma_tuning = 1e-6)


fit_gddlm05 <- gddlm(mu_formula = make_mu_formula,
                     sigma_formula = make_sigma_formula,
                     data = data.frame(response = Y, X[,-1]),
                     gamma_tuning = 0.5)


fit_RobGam <- DoubleRobGam(formulaM = make_mu_formula,
                           formulaG = make_sigma_formula,
                           selection = "none",
                           data = data.frame(response = Y, X[,-1]))


fit_psidiv <- try(GJRM::gamlss(list(make_mu_formula, make_sigma_formula),
                               data = data.frame(response = Y, X[,-1]),
                               family = "N",
                               robust = TRUE),
                  silent = TRUE) #' This takes the longest by far!


data.frame(nonrobust = fit_lmls$coefficients[1:(p+1)],
           gddlm_opt = fit_gddlm_choosegamma$mu_coefficients,
           gddlm_gamma0 = fit_gddlm0$mu_coefficients,
           gddlm_gamm05 = fit_gddlm05$mu_coefficients,
           Huber = fit_RobGam$fitM$coefficients,
           Psidiv = fit_psidiv$coefficients[1:(p+1)])

data.frame(nonrobust = fit_lmls$coefficients[-(1:(p+1))],
           gddlm_opt = fit_gddlm_choosegamma$sigma_coefficients,
           gddlm_gamma0 = fit_gddlm0$sigma_coefficients,
           gddlm_gamm05 = fit_gddlm05$sigma_coefficients,
           Huber = fit_RobGam$fitG$coefficients*0.5,
           Psidiv = fit_psidiv$coefficients[-(1:(p+1))])



# Attempt variable selection
fit_gddlm_choosegama_lsa <- gddlm_lsa(object = fit_gddlm_choosegamma)
fit_gddlm_gamma0_lsa <- gddlm_lsa(object = fit_gddlm0)
fit_gddlm_gamma05_lsa <- gddlm_lsa(object = fit_gddlm05)
dat <- data.frame(y = Y, X[,-1])
fit_gamlsslasso <- gamlss(y ~ gnet(x.vars = names(dat)[-1]),
                          sigma.fo = ~ gnet(x.vars = names(dat)[-1]),
                          data = dat,
                          family = NO,
                          i.control = glim.control(cyc=1, bf.cyc=1))
fit_gamlsslasso



robgam_coef <- fit_RobGam$fitM$coefficients
robgam_coef[abs(fit_RobGam$fitM$coefficients / sqrt(diag(fit_RobGam$fitM$cov.coef))) < qnorm(0.975)] <- 0
psidiv_coef <- fit_psidiv$coefficients[1:(p+1)]
psidiv_coef[abs(fit_psidiv$coefficients[1:(p+1)] / sqrt(diag(fit_psidiv$Vb))[1:(p+1)]) < qnorm(0.975)] <- 0
data.frame(
    true = Beta_true,
    gddlm_opt = fit_gddlm_choosegama_lsa$mu_coefficients_path[,which.min(fit_gddlm_choosegama_lsa$information$BIC)],
    gddlm_0 = fit_gddlm_gamma0_lsa$mu_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)],
    gddlm_05 = fit_gddlm_gamma05_lsa$mu_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)],
    Huber = robgam_coef,
    psidiv = psidiv_coef,
    gamlss_lasso = c(NA, tail(getSmo(fit_gamlsslasso, "mu") ,1)[[1]]$beta))
rm(robgam_coef, psidiv_coef)



robgam_coef <- fit_RobGam$fitG$coefficients*0.5
psidiv_coef <- fit_psidiv$coefficients[-(1:(p+1))]
psidiv_coef[abs(fit_psidiv$coefficients[-(1:(p+1))] / sqrt(diag(fit_psidiv$Vb))[-(1:(p+1))])] <- 0
data.frame(
    true = Delta_true,
    gddlm_opt = fit_gddlm_choosegama_lsa$sigma_coefficients_path[,which.min(fit_gddlm_choosegama_lsa$information$BIC)],
    gddlm_0 = fit_gddlm_gamma0_lsa$sigma_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)],
    gddlm_05 = fit_gddlm_gamma05_lsa$sigma_coefficients_path[,which.min(fit_gddlm_gamma0_lsa$information$BIC)],
    Huber = robgam_coef,
    psidiv = psidiv_coef,
    gamlss_lasso = c(NA, tail(getSmo(fit_gamlsslasso, "sigma") ,1)[[1]]$beta))
rm(robgam_coef, psidiv_coef)




#----------------------
# Custom testing
#----------------------
function() {
    mu_formula = make_mu_formula
    sigma_formula = make_sigma_formula
    data = data.frame(response = Y, X[,-1])
    gamma_tuning = 1e-12
    mu_bar_lambda = 0.2
    sigma_bar_lambda = 0.2
    control = list(max_iteration = 1e3, tol = 1e-5, dampening = c(1,1), trace = FALSE)


    object = fit_gddlm_choosegamma
    nlambda = 100
    lambda_min_ratio = 1e-6
    lambda = NULL
    min_df = 0
    num_cores = NULL
    prior_weights = NULL
    control = list(maxit = 1000, eps = 1e-6, round_eps = 1e-6)
    }
