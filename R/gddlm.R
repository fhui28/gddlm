#' @title Fast robust double linear models using gamma divergence.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function fits a robust double linear model, where both the mean (mu) and the log of the standard deviation (sigma) are modeled as a function of predictors, based on applying the [gamma divergence](https://www.mdpi.com/1099-4300/19/11/608) to the associated normally distributed likelihood function (hence the name gddlm, which stands for gamma divergence double linear models).
#'
#' Optimization is performed via a minorization-maximization (MM) algorithm, while standard errors for the coefficient estimates in the mean and standard deviation are constructed based on the inverse of the observed Fisher information matrix for the gamma divergence function.
#'
#'
#' @param mu_formula A double-sided formula describing how the response is regressed against the covariates in the mean part of the double linear model.
#' @param sigma_formula A single-sided formula describing how the log of the standard deviation is regressed against the covariates in the double linear model.
#' @param data A data frame containing the response and covariates used in `mu_formula` and `sigma_formula`.
#' @param gamma_tuning A positive tuning parameter for the gamma divergence, and typically less than 1 for fitting robust double linear models. The closer this value is to zero, the more the gamma divergence behaves like standard maximum likelihood estimation e.g., the same as what is done in [mgcv::gam()] with the [mgcv::gaulss()] family, or in the [lmls::lmls()] function.
#' @param control A list of control parameters for fitting the gamma divergence double linear model. The list can contain the following elements:
#' \itemize{
#' \item \code{max_iteration}: The maximum number of iterations for the MM algorithm.
#' \item \code{tol}: The tolerance for convergence of the MM algorithm.
#' \item \code{dampening}: A two element vector that can be used to dampen the starting values of the coefficients in the mean and standard deviation parts of the double linear model. By default, the function obtains starting values from a non-robust double linear model fit, and it is possible these might be poor starting values e.g., the magnitudes of the coefficients might be too large and the MM algorithm subsequently struggles. To counter this, dampening values less than one can be used to "dampen" or reduce the magnitudes of these starting values, which may help to improve convergence. The first element in the vector controls the mean part, while the second element controls the standard deviation part of the double linear model. Defaults to c(1, 1), which corresponds to no dampening.
#' \item \code{trace} Should the function print out trace information during the MM algorithm?
#' }
#'
#'
#' @details
#' Please see manuscript for more details (currently in preparation).
#'
#' @return An object of class `gddlm` with the following components:
#' \itemize{
#' \item \code{call} The matched call.
#' \item \code{mu_formula} As per input argument.
#' \item \code{sigma_formula} As per input argument.
#' \item \code{model} Same as `data` input argument.
#' \item \code{mu_coefficients} The estimated coefficients for the mean part of the double linear model.
#' \item \code{sigma_coefficients} The estimated coefficients for the log standard deviation part of the double linear model.
#' \item \code{mu_fitted} The fitted values for the mean in the double linear model.
#' \item \code{sigma_fitted} The fitted values for the standard deviation in the double linear model.
#' \item \code{raw_residuals} The raw residuals from the double linear model, which is the difference between the response and the fitted values for the mean.
#' #' \item \code{pearson_residuals} The Pearson residuals from the double linear model, which is the raw residuals divided by the fitted values for the standard deviation.
#' \item \code{final_weights} The final weights from the MM algorithm. The weights sum to the total number of observations in the data. Weights comparably less than 1 indicate that the observation is down-weighted in the estimation process i.e., is considered more likely to be an outlying observation in the double linear model.
#' \item \code{gamma_tuning} As per input argument.
#' \item \code{gamma_divergence} The value of the gamma divergence at convergence.
#' \item \code{V_mu_coefficients} The estimated covariance (inverse observed information) matrix for the coefficients corresponding to the mean part of the double linear model.
#' \item \code{V_sigma_coefficients} The estimated covariance (inverse observed information) matrix for the coefficients corresponding to the log standard deviation part of the double linear model.
#' \item \code{V_full} The full estimated covariance (inverse observed information) matrix.
#' }
#'
#'
#' @author
#' Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(GJRM)
#' library(mgcv)
#' library(mvtnorm)
#'
#' # Simulation setting and data generation
#' set.seed(012026)
#' n <- 250
#' p <- 10
#' H <- abs(outer(1:p, 1:p, "-"))
#' H <- 0.75^H
#' X <- cbind(1, rmvnorm(n, sigma = H))
#' Beta_true <- c(-1, -1, 1, -0.5, 0.5, 0.2, -0.2, rep(0, p - 6))
#' Delta_true <- c(0.5, 0.5, 0, 0, -0.5, 0, 0, 0.25, -0.25, rep(0, p - 8))
#' ## Generate responses
#' mu_true <- as.vector(X %*% Beta_true)
#' sig_true <- exp(as.vector(X %*% Delta_true))
#' Y <- rnorm(n, mean = mu_true, sd = sig_true)
#'
#' ## Create outliers
#' om <- 0.05
#' out <- rbinom(n, 1, om)
#' Y[out==1] <- Y[out==1] + 50
#'
#' hist(Y, breaks = 20)
#'
#'
#' # Fit gamma divergence double linear models
#' make_mu_formula <- paste0("response ~" , paste0("X", 1:p, collapse = " + ")) %>%
#' formula
#' make_sigma_formula <- paste0("~" , paste0("X", 1:p, collapse = " + ")) %>%
#' formula
#'
#' ## Select the gamma tuning parameter via cross-validation
#' runcv <- cv_gddlm(mu_formula = make_mu_formula,
#' sigma_formula = make_sigma_formula,
#' data = data.frame(response = Y, X[,-1]))
#'
#' fit_gddglm_choosegamma <- gddlm(mu_formula = make_mu_formula,
#' sigma_formula = make_sigma_formula,
#' data = data.frame(response = Y, X[,-1]),
#' gamma_tuning = runcv$gamma_optimal)
#'
#' # Demonstration of some of the functions available in the package
#' fit_gddglm_choosegamma
#'
#' summary(fit_gddglm_choosegamma)
#'
#' gratia::appraise(fit_gddglm_choosegamma)
#'
#' coef.gddlm(fit_gddglm_choosegamma)
#' coef.gddlm(fit_gddglm_choosegamma, type = "sigma")
#'
#' fitted.gddlm(fit_gddglm_choosegamma, type = "sigma")
#' residuals.gddlm(fit_gddglm_choosegamma, type = "pearson") %>%
#' hist(breaks = 20, main = "Pearson residuals")
#'
#' vcov.gddlm(fit_gddglm_choosegamma)
#' vcov.gddlm(fit_gddglm_choosegamma, type = "sigma")
#'
#' simulate(fit_gddglm_choosegamma, nsim = 10)
#'
#' predict(fit_gddglm_choosegamma, predictor = "mu", se.fit = TRUE)
#' predict(fit_gddglm_choosegamma, predictor = "sigma", se.fit = TRUE)
#'
#'
#' # Fitting some other double linear models
#' fit_gddlmgamma01 <- gddlm(mu_formula = make_mu_formula,
#' sigma_formula = make_sigma_formula,
#' data = data.frame(response = Y, X[,-1]),
#' gamma_tuning = 0.1)
#'
#' fit_gddlmgamma02 <- gddlm(mu_formula = make_mu_formula,
#' sigma_formula = make_sigma_formula,
#' data = data.frame(response = Y, X[,-1]),
#' gamma_tuning = 0.2)
#'
#' fit_gddlmgamma05 <- gddlm(mu_formula = make_mu_formula,
#' sigma_formula = make_sigma_formula,
#' data = data.frame(response = Y, X[,-1]),
#' gamma_tuning = 0.5)
#'
#' ## Non-robust approach
#' fit_lmls <- gam(list(make_mu_formula, make_sigma_formula),
#' data = data.frame(response = Y, X[,-1]),
#' family = gaulss(),
#' method = "ML")
#'
#' ## Psi-divergence. Note this takes the longest by far...
#' fit_psidiv <- GJRM::gamlss(list(make_mu_formula, make_sigma_formula),
#' data = data.frame(response = Y, X[,-1]),
#' family = "N",
#' robust = TRUE)
#'
#'
#' ## Assess performance
#' Beta_est <- data.frame(truevalue = Beta_true,
#' lmls = fit_lmls$coefficients[1:ncol(X)],
#' gamma01 = fit_gddlmgamma01$mu_coefficients,
#' gamma02 = fit_gddlmgamma02$mu_coefficients,
#' gamma05 = fit_gddlmgamma05$mu_coefficients,
#' choosegamma = fit_gddglm_choosegamma$mu_coefficients,
#' psidiv = fit_psidiv$coefficients[1:ncol(X)])
#'
#' Delta_est <- data.frame(truevalue = Delta_true,
#' lmls = fit_lmls$coefficients[-(1:ncol(X))],
#' gamma01 = fit_gddlmgamma01$sigma_coefficients,
#' gamma02 = fit_gddlmgamma02$sigma_coefficients,
#' gamma05 = fit_gddlmgamma05$sigma_coefficients,
#' choosegamma = fit_gddglm_choosegamma$sigma_coefficients,
#' psidiv = fit_psidiv$coefficients[-(1:ncol(X))])
#'
#'
#' Beta_est
#' Beta_est %>% #' Bias averaged across coefficients (note only one dataset!)
#' mutate(across(-truevalue, ~ .x - truevalue)) %>%
#' select(-truevalue) %>%
#' colMeans
#'
#' Beta_est %>%
#' mutate(across(-truevalue, ~ (.x - truevalue)^2)) %>%
#' select(-truevalue) %>%
#' colSums #' Frobenius norm between estimate and true coefficients (not only one dataset!)
#'
#' Delta_est
#' Delta_est %>% #' Bias averaged across coefficients (not only one dataset!)
#' mutate(across(-truevalue, ~ .x - truevalue)) %>%
#' select(-truevalue) %>%
#' colMeans
#'
#' Delta_est %>%
#' mutate(across(-truevalue, ~ (.x - truevalue)^2)) %>%
#' select(-truevalue) %>%
#' colSums #' Frobenius norm between estimate and true coefficients (not only one dataset!)
#' }
#'
#' @export gddlm
#' @importFrom Matrix Diagonal
#' @importFrom mgcv gam model.matrix.gam gaulss
#' @importFrom stats dnorm gaussian model.frame model.response update.formula
#' @md
gddlm <- function(mu_formula,
                  sigma_formula,
                  data,
                  gamma_tuning = 0.1,
                  control = list(max_iteration = 1e3, tol = 1e-5, dampening = c(1,1), trace = FALSE)) {

    ##-----------------------
    #' # Checks and balances
    ##-----------------------
    if(is.null(control$max_iteration))
        control$max_iteration <- 1e3
    if(is.null(control$tol))
        control$tol <- 1e-4
    if(is.null(control$dampening))
        control$dampening <- c(1,1)
    if(is.null(control$trace))
        control$trace <- FALSE

    n <- nrow(data)
    fake_gam_mu <- mgcv::gam(mu_formula, data = data, family = gaussian())
    mu_X <- mgcv::model.matrix.gam(fake_gam_mu)
    Y <- model.frame(fake_gam_mu)
    Y <- as.vector(model.response(Y))
    rm(fake_gam_mu)
    expand_sigma_formula <- update.formula(mu_formula, sigma_formula)
    fake_gam_sig <- mgcv::gam(expand_sigma_formula, data = data, family = gaussian())
    sig_X <- mgcv::model.matrix.gam(fake_gam_sig)
    rm(fake_gam_sig, expand_sigma_formula)

    fake_gam_lmls <- mgcv::gam(list(mu_formula, sigma_formula),
                               data = data,
                               family = gaulss())

    ##-----------------------
    #' # initial values
    ##-----------------------
    fit_start <- gam(list(mu_formula, sigma_formula),
                     data = data,
                     family = gaulss(), method = "ML")
    mu_coefficients <- fit_start$coefficients[1:ncol(mu_X)]^control$dampening[1]
    sigma_coefficients <- fit_start$coefficients[-(1:ncol(mu_X))]^control$dampening[2]
    old_gamma_divergence <- -Inf


    ##-----------------------
    #' # Run MM-algorithm
    ##-----------------------
    for(r in 1:control$max_iteration) {
        # Minorization -- construct density weights
        mu_coefficients_old <- mu_coefficients
        sigma_coefficients_old <- sigma_coefficients
        mu <- as.vector(mu_X %*% mu_coefficients)
        sig <- exp(as.vector(sig_X %*% sigma_coefficients))
        dens <- dnorm(Y, mean = mu, sd = sig)
        ww <- dens^gamma_tuning / mean(dens^gamma_tuning)


        # UPDATING METHOD VERSION 2 -- *CURRENTLY DOES NOT WORK AS GAULSS FITTING DOES NOT PERMIT WEIGHTS*
        # fit_update <- gam(list(mu_formula, sigma_formula),
        #                   data = data,
        #                   weights = ww,
        #                   family = gaulss(),
        #                   method = "ML")
        # mu_coefficients <- fit_update$coefficients[1:ncol(mu_X)]
        # sigma_coefficients <- fit_update$coefficients[-(1:ncol(mu_X))]


        # UPDATING METHOD VERSION 1
        # update mu_coefficients
        H <- crossprod(mu_X, ww * mu_X / sig^2)
        mu_coefficients <- solve(H) %*% crossprod(mu_X, ww * Y / sig^2)

        # update sigma_coefficients -- Using Fisher scoring here instead of straight up Newton-Raphson
        resid <- Y - as.vector(mu_X %*% mu_coefficients)
        z <- mean(1/sig^gamma_tuning)
        sig_weights <- 2 * ww + gamma_tuning^2/((1 + gamma_tuning) * z) * 1/sig^gamma_tuning
        H <- crossprod(sig_X, sig_weights * sig_X) # Expected Fisher information
        Grad <- apply((ww * resid^2 / sig^2 - ww) * sig_X + (gamma_tuning^2/((1 + gamma_tuning) * z) * 1/sig^gamma_tuning) * sig_X , 2, sum)
        sigma_coefficients <- sigma_coefficients + solve(H) %*% Grad

        # Convergence check
        new_gamma_div <- n / gamma_tuning * log(mean(dens^gamma_tuning)) - n / (1 + gamma_tuning) * log(mean(1/sig^gamma_tuning))
        new_gamma_div <- new_gamma_div
        dd <- new_gamma_div - old_gamma_divergence
        if(r %% 10 == 0 & control$trace == TRUE)
            message("Iteration: ", r, "\t New gamma divergence: ", round(new_gamma_div, 4),
                    "\t Old gamma divergence: ", round(old_gamma_divergence, 4),
                    "\t Difference: ", round(dd, 4))

        old_gamma_divergence <- new_gamma_div
        if(dd < control$tol)
            break
        }


    ##-----------------------
    #' # Hessian (mu_coefficients and sigma_coefficients)
    ##-----------------------
    mu <- as.vector(mu_X %*% mu_coefficients)
    sig <- exp(as.vector(sig_X %*% sigma_coefficients))
    ww <- dnorm(Y, mean = mu, sd = sig)^gamma_tuning
    S <- mean(ww)
    Z <- (Y - mu) / sig^2
    A <- apply(mu_X * ww * Z, 2, sum)
    a <- ww * (gamma_tuning * Z^2 - 1 / sig^2)
    Hess_mu_coefficients <- (-1) * gamma_tuning / (n * S^2) * tcrossprod(A) + (1/S) * crossprod(mu_X, a * mu_X)

    tau <- (Y - mu)^2 / sig^2 - 1
    a <- ww * (gamma_tuning * tau^2 - 2 * (tau + 1))
    A <- apply(sig_X * ww * tau, 2, sum)
    c_i <- 1/sig^gamma_tuning / sum(1/sig^gamma_tuning)
    Sigma2 <- crossprod(sig_X, c_i * sig_X) - tcrossprod(apply(c_i * sig_X, 2, sum))
    Hess_sigma_coefficients <- (-1) * gamma_tuning / (n * S^2) * tcrossprod(A) + (1 / S) * crossprod(sig_X, a * sig_X) -  n * gamma_tuning^2/(1 + gamma_tuning) * Sigma2

    Hess_musigma_coefficients <- (1 / S) * crossprod(mu_X, ww * (gamma_tuning * tau * Z - 2 * Z) * sig_X) - gamma_tuning / (n * S^2) * tcrossprod(apply(mu_X * ww * Z, 2, sum), apply(sig_X * ww * tau, 2, sum))


    ##-----------------------
    #' # Final wrap up
    ##-----------------------
    # C <- ( (1 + gamma_tuning)^(-1/2) * (2*pi*sig^2)^(-gamma_tuning/2) )^(gamma_tuning / (1 + gamma_tuning))
    # dens <- dnorm(Y, mean = mu, sd = sig)
    # term1 <- dens^(gamma_tuning) * 2*(gamma_tuning * (Y - mu)^2 - sig^2) / ( sig^4 * C)
    # term2 <- dens^(2*gamma_tuning) * (Y - mu)^2 / ( sig^4 * C^2)
    # Hscore <- mean(term1 + term2)

    dens <- dnorm(Y, mean = mu, sd = sig)
    ww <- dens^gamma_tuning / mean(dens^gamma_tuning)

    fullinfo <- -rbind(cbind(Hess_mu_coefficients, Hess_musigma_coefficients),
                      cbind(t(Hess_musigma_coefficients), Hess_sigma_coefficients))

    out <- list(call = match.call(),
                mu_formula = mu_formula,
                sigma_formula = sigma_formula,
                model = data,
                mu_coefficients = as.vector(mu_coefficients),
                sigma_coefficients = as.vector(sigma_coefficients),
                mu_fitted = mu,
                sigma_fitted = sig,
                raw_residuals = Y - mu,
                pearson_residuals = (Y - mu) / sig,
                final_weights = ww,
                gamma_tuning = gamma_tuning,
                gamma_divergence = new_gamma_div,
                V_full = solve(fullinfo))
    out$V_mu_coefficients <- out$V_full[1:ncol(mu_X), 1:ncol(mu_X)]
    out$V_sigma_coefficients <- out$V_full[-(1:ncol(mu_X)), -(1:ncol(mu_X))]
    names(out$mu_coefficients) <- colnames(mu_X)
    names(out$sigma_coefficients) <- colnames(sig_X)
    rownames(out$V_mu_coefficients) <- colnames(out$V_mu_coefficients) <- colnames(mu_X)
    rownames(out$V_sigma_coefficients) <- colnames(out$V_sigma_coefficients) <- colnames(sig_X)
    rownames(out$V_full) <- colnames(out$V_full) <- c(colnames(mu_X), colnames(sig_X))

    class(out) <- "gddlm"
    return(out)
    }



#'------------------
#' **OLD CODE THAT ATTEMPTED TO INCORPORATE PENALIZATION VIA THE BAR PENALTY**
#' Abandoned as it is did not work very well and tuning parameter became too much of an issue! It seems better to separate it out into separate stages, and do an robust LQA?
#'------------------

#' gddlm <- function(mu_formula,
#'                   sigma_formula,
#'                   data,
#'                   gamma_tuning = 0.1,
#'                   mu_bar_lambda = 0,
#'                   sigma_bar_lambda = 0,
#'                   control = list(max_iteration = 1e3, tol = 1e-5, dampening = c(1,1), trace = FALSE),
#'                   bar_control = list(tol = 1e-4, max_iteration = 1, round_eps = 1e-6)) {
#'
#'     ##-----------------------
#'     #' # Checks and balances
#'     ##-----------------------
#'     if(is.null(control$max_iteration))
#'         control$max_iteration <- 1e3
#'     if(is.null(control$tol))
#'         control$tol <- 1e-4
#'     if(is.null(control$dampening))
#'         control$dampening <- c(1,1)
#'     if(is.null(control$trace))
#'         control$trace <- FALSE
#'
#'     if(is.null(bar_control$tol))
#'         bar_control$tol <- 1e-4
#'     if(is.null(bar_control$max_iteration))
#'         bar_control$max_iteration <- 1e3
#'     if(is.null(bar_control$round_eps))
#'         bar_control$round_eps <- 1e-6
#'
#'
#'     n <- nrow(data)
#'     fake_gam_mu <- mgcv::gam(mu_formula, data = data, family = gaussian())
#'     mu_X <- mgcv::model.matrix.gam(fake_gam_mu)
#'     Y <- model.frame(fake_gam_mu)
#'     Y <- as.vector(model.response(Y))
#'     rm(fake_gam_mu)
#'     expand_sigma_formula <- update.formula(mu_formula, sigma_formula)
#'     fake_gam_sig <- mgcv::gam(expand_sigma_formula, data = data, family = gaussian())
#'     sig_X <- mgcv::model.matrix.gam(fake_gam_sig)
#'     rm(fake_gam_sig, expand_sigma_formula)
#'
#'     fake_gam_lmls <- mgcv::gam(list(mu_formula, sigma_formula),
#'                                data = data,
#'                                family = gaulss())
#'
#'     ##-----------------------
#'     #' # initial values
#'     ##-----------------------
#'     fit_start <- gam(list(mu_formula, sigma_formula),
#'                      data = data,
#'                      family = gaulss(), method = "ML")
#'     mu_coefficients <- fit_start$coefficients[1:ncol(mu_X)]^control$dampening[1]
#'     sigma_coefficients <- fit_start$coefficients[-(1:ncol(mu_X))]^control$dampening[2]
#'     old_gamma_divergence <- -Inf
#'
#'
#'     ##-----------------------
#'     #' # Run MM-algorithm
#'     ##-----------------------
#'     for(r in 1:control$max_iteration) {
#'         # Minorization -- construct density weights
#'         mu_coefficients_old <- mu_coefficients
#'         sigma_coefficients_old <- sigma_coefficients
#'         mu <- as.vector(mu_X %*% mu_coefficients)
#'         sig <- exp(as.vector(sig_X %*% sigma_coefficients))
#'         dens <- dnorm(Y, mean = mu, sd = sig)
#'         ww <- dens^gamma_tuning / mean(dens^gamma_tuning)
#'
#'
#'         # UPDATING METHOD VERSION 2 -- *CURRENTLY DOES NOT WORK AS GAULSS FITTING IS NOT WEIGHTED*
#'         # fit_update <- gam(list(mu_formula, sigma_formula),
#'         #                   data = data,
#'         #                   weights = ww,
#'         #                   family = gaulss(),
#'         #                   method = "ML")
#'         # mu_coefficients <- fit_update$coefficients[1:ncol(mu_X)]
#'         # sigma_coefficients <- fit_update$coefficients[-(1:ncol(mu_X))]
#'
#'
#'         # UPDATING METHOD VERSION 1
#'         # update mu_coefficients
#'         err <- Inf
#'         counter <- 0
#'         cw_mu_coefficients <- mu_coefficients
#'         Dbar <- Diagonal(n = length(cw_mu_coefficients))
#'         Dbar[1,1] <- 0 #' Do not penalize intercept
#'         while(counter < bar_control$max_iteration & err > bar_control$tol) {
#'             GammaMatrix <- Diagonal(x = cw_mu_coefficients, n = length(cw_mu_coefficients))
#'             GammaMatrix[1,1] <- 1 #' Do not penalize intercept
#'
#'             H <- crossprod(mu_X, ww * mu_X / sig^2)
#'             mu_coefficients <- GammaMatrix %*% solve(GammaMatrix %*% H %*% GammaMatrix + mu_bar_lambda * Dbar) %*% GammaMatrix %*% crossprod(mu_X, ww * Y / sig^2)
#'             mu_coefficients <- as.vector(mu_coefficients)
#'
#'             err <- sum((mu_coefficients - cw_mu_coefficients)^2)
#'             cw_mu_coefficients <- mu_coefficients
#'             counter <- counter + 1
#'         }
#'
#'         #' Note the convergence properties of the algorithm does basically ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested. In this case, I will round based on eps = 1e-6
#'         mu_coefficients[abs(mu_coefficients) < bar_control$round_eps] <- 0
#'         rm(err, counter, cw_mu_coefficients, Dbar, GammaMatrix, H)
#'
#'
#'         # update sigma_coefficients -- Using Fisher scoring here instead of straight up Newton-Raphson
#'         err <- Inf
#'         counter <- 0
#'         resid <- Y - as.vector(mu_X %*% mu_coefficients)
#'         cw_sigma_coefficients <- sigma_coefficients
#'         Dbar <- Diagonal(n = length(cw_sigma_coefficients))
#'         Dbar[1,1] <- 0 #' Do not penalize intercept
#'         Dvec <- matrix(1, nrow = length(cw_sigma_coefficients), ncol = 1)
#'         Dvec[1,1] <- 0 #' Do not penalize intercept
#'         while(counter < bar_control$max_iteration & err > bar_control$tol) {
#'             GammaMatrix <- Diagonal(x = cw_sigma_coefficients, n = length(cw_sigma_coefficients))
#'             GammaMatrix[1,1] <- 1 #' Do not penalize intercept
#'
#'             H <- 2 * crossprod(sig_X, ww * sig_X) # Expected Fisher information
#'             #Grad <- apply(ww * resid^2 * sig_X / sig^2 - sig_X / (1 + gamma_tuning), 2, sum) - sigma_bar_lambda * W %*% cw_sigma_coefficients
#'             GammaMatrixGrad <- GammaMatrix %*% apply(ww * resid^2 * sig_X / sig^2 - sig_X / (1 + gamma_tuning), 2, sum) - sigma_bar_lambda * Dvec
#'
#'             sigma_coefficients <- cw_sigma_coefficients + GammaMatrix %*% solve(GammaMatrix %*% H %*% GammaMatrix + sigma_bar_lambda * Dbar) %*% GammaMatrixGrad
#'             sigma_coefficients <- as.vector(sigma_coefficients)
#'
#'             err <- sum((sigma_coefficients - cw_sigma_coefficients)^2)
#'             cw_sigma_coefficients <- sigma_coefficients
#'             sig <- exp(as.vector(sig_X %*% cw_sigma_coefficients))
#'             counter <- counter + 1
#'         }
#'
#'         #' Note the convergence properties of the algorithm does basically ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested. In this case, I will round based on eps = 1e-6
#'         sigma_coefficients[abs(sigma_coefficients) < bar_control$round_eps] <- 0
#'         rm(err, counter, cw_sigma_coefficients, Dbar, GammaMatrix, H, Dvec, resid, GammaMatrixGrad)
#'
#'
#'         # Convergence check
#'         new_gamma_div <- n / gamma_tuning * log(mean(dens^gamma_tuning)) + gamma_tuning / (2 * (1 + gamma_tuning)) * sum(log(sig^2))
#'         #make_quad_penalty <- 0.5 * as.vector(t(mu_coefficients) %*% S_mu %*% mu_coefficients) + 0.5 * as.vector(t(sigma_coefficients) %*% S_sigma %*% sigma_coefficients)
#'         new_gamma_div <- new_gamma_div #- make_quad_penalty
#'         dd <- new_gamma_div - old_gamma_divergence
#'         if(r %% 10 == 0 & control$trace == TRUE)
#'             message("Iteration: ", r, "\t New gamma divergence: ", round(new_gamma_div, 4),
#'                     "\t Old gamma divergence: ", round(old_gamma_divergence, 4),
#'                     "\t Difference: ", round(dd, 4))
#'
#'         old_gamma_divergence <- new_gamma_div
#'         if(dd < control$tol)
#'             break
#'     }
#'
#'
#'     ##-----------------------
#'     #' # Hessian (mu_coefficients and sigma_coefficients)
#'     ##-----------------------
#'     mu <- as.vector(mu_X %*% mu_coefficients)
#'     sig <- exp(as.vector(sig_X %*% sigma_coefficients))
#'     ww <- dnorm(Y, mean = mu, sd = sig)^gamma_tuning
#'     S <- mean(ww)
#'     Z <- (Y - mu) / sig^2
#'     A <- apply(mu_X * ww * Z, 2, sum)
#'     a <- ww * (gamma_tuning * Z^2 - 1 / sig^2)
#'     Hess_mu_coefficients <- (-1) * gamma_tuning / (n * S^2) * tcrossprod(A) + (1/S) * crossprod(mu_X, a * mu_X)
#'
#'     Z <- (Y - mu)^2 / sig^2 - 1
#'     a <- ww * (0.5 * gamma_tuning * Z^2 - (Z + 1))
#'     A <- apply(sig_X * ww * Z, 2, sum)
#'     Hess_sigma_coefficients <- (-1) * gamma_tuning / (n * S^2) * tcrossprod(A) + (2 / S) * crossprod(sig_X, a * sig_X)
#'
#'
#'     ##-----------------------
#'     #' # Final wrap up
#'     ##-----------------------
#'     # C <- ( (1 + gamma_tuning)^(-1/2) * (2*pi*sig^2)^(-gamma_tuning/2) )^(gamma_tuning / (1 + gamma_tuning))
#'     # dens <- dnorm(Y, mean = mu, sd = sig)
#'     # term1 <- dens^(gamma_tuning) * 2*(gamma_tuning * (Y - mu)^2 - sig^2) / ( sig^4 * C)
#'     # term2 <- dens^(2*gamma_tuning) * (Y - mu)^2 / ( sig^4 * C^2)
#'     # Hscore <- mean(term1 + term2)
#'
#'     dens <- dnorm(Y, mean = mu, sd = sig)
#'     ww <- dens^gamma_tuning / mean(dens^gamma_tuning)
#'
#'     out <- list(call = match.call(),
#'                 mu_formula = mu_formula,
#'                 sigma_formula = sigma_formula,
#'                 mu_bar_lambda = mu_bar_lambda,
#'                 sigma_bar_lambda = sigma_bar_lambda,
#'                 model = fake_gam_lmls$model,
#'                 mu_coefficients = as.vector(mu_coefficients),
#'                 sigma_coefficients = as.vector(sigma_coefficients),
#'                 mu_fitted = mu,
#'                 sigma_fitted = sig,
#'                 raw_residuals = Y - mu,
#'                 pearson_residuals = (Y - mu) / sig,
#'                 final_weights = ww,
#'                 gamma_tuning = gamma_tuning,
#'                 gamma_divergence = new_gamma_div,
#'                 V_mu_coefficients = solve(-Hess_mu_coefficients),
#'                 V_sigma_coefficients = solve(-Hess_sigma_coefficients))
#'     names(out$mu_coefficients) <- colnames(mu_X)
#'     names(out$sigma_coefficients) <- colnames(sig_X)
#'     rownames(out$V_mu_coefficients) <- colnames(out$V_mu_coefficients) <- colnames(mu_X)
#'     rownames(out$V_sigma_coefficients) <- colnames(out$V_sigma_coefficients) <- colnames(sig_X)
#'
#'     class(out) <- "gddlm"
#'     return(out)
#' }
