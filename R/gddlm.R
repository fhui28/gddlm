#' @title Fast robust double linear models using gamma divergence.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function fits a robust double linear model, where both the mean (mu) and the log of the standard deviation (sigma) are modeled as a function of predictors, based on applying the gamma divergence to the associated normally distributed likelihood function (hence the name gddlm, which stands for gamma divergence double linear models).
#'
#' Optimization is performed via a minorization-maximization (MM) algorithm, while standard errors for the coefficient estimates in the mean and standard deviation are constructed based on the inverse of the observed Fisher information matrix for the gamma divergence function.
#'
#'
#' @param mu_formula A double-sided formula describing how the response is regressed against the covariates in the mean part of the double linear model. 
#' @param sigma_formula A single-sided formula describing how the log of the standard deviation is regressed against the covariates in the double linear model. 
#' @param data A data frame containing the response and covariates used in \code{mu_formula} and \code{sigma_formula}.
#' @param gamma_tuning A positive tuning parameter for the gamma divergence, and typically less than 1 for fitting robust double linear models. The closer this value is to zero, the more the gamma divergence behaves like standard maximum likelihood estimation e.g., the same as what is done in [mgcv::gam()] with the [mgcv::gaulss()] family, or in the [lmls::lmls()] function.
#' @param max_iteration The maximum number of iterations for the MM algorithm. Default is 1000.
#' @param tol The tolerance for convergence of the MM algorithm. Default is 1e-6.
#' @param dampening A two element vector that can be used to dampen the starting values of the coefficients in the mean and standard deviation parts of the double linear model. By default, the function obtains starting values from a non-robust double linear model fit, and it is possible these might be poor starting values e.g., the magnitudes of the coefficients might be too large and the MM algorithm subsequently struggles. To counter this, dampening values less than one can be used to "dampen" or reduce the magnitudes of these starting values, which may help to improve convergence. The first element in the vector controls the mean part, while the second element controls the standard deviation part of the double linear model. Defaults to c(1, 1), which corresponds to no dampening.
#' 
#' 
#' @details
#' See manuscript for more details (currently in preparation).
#'
#' @return An object of class `gddlm` with the following components:
#' \itemize{
#' \item \code{call} The matched call.
#' \item \code{mu_coefficients} The estimated coefficients for the mean part of the double linear model.
#' \item \code{sigma_coefficients} The estimated coefficients for the log standard deviation part of the double linear model.
#' \item \code{mu_fitted} The fitted values for the mean in the double linear model. 
#' \item \code{sigma_fitted} The fitted values for the standard deviation in the double linear model. 
#' \item \code{raw_residuals} The raw residuals from the double linear model, which is the difference between the response and the fitted values for the mean.
#' #' \item \code{pearson_residuals} The Pearson residuals from the double linear model, which is the raw residuals divided by the fitted values for the standard deviation.
#' \item \code{final_weights} The final weights from the MM algorithm. The weights sum to the total number of observations in the data. Weights comparably less than 1 indicate that the observation is down-weighted in the estimation process i.e., is considered more likely to be an outlying observation in the double linear model.
#' \item \code{gamma_tuning} The gamma tuning parameter used in the estimation process.
#' \item \code{gamma_divergence} The value of the gamma divergence at convergence.
#' \item \code{Hscore} The Hyvarinen score (H-score) of the double linear model. This is a measure which can be used to select the tuning parameter `gamma_tuning`, when comparing across difference `gdam` fits. The H-score calculation is based on the proposal of [Sugusawa and Yonekura, 2021](https://doi.org/10.3390/e23091147).
#' \item \code{V_mu_coefficients} The estimated covariance (inverse Fisher information) matrix for the coefficients corresponding to the mean part of the double linear model.
#' \item \code{V_sigma_coefficients} The estimated covariance (inverse Fisher information) matrix for the coefficients corresponding to the log standard deviation part of the double linear model.
#' }
#' 
#' 
#' @author
#' Francis K.C. Hui <francis.hui@anu.edu.au>
#' 
#' @examples
#' \dontrun{
#' TO BE DONE!!!
#' }
#' @export gddlm
#' @importFrom mgcv gam model.matrix.gam gaulss
#' @importFrom stats dnorm 
#' @md
gddlm <- function(mu_formula,
                  sigma_formula,
                  data,
                  gamma_tuning = 0.1, 
                  max_iteration = 1e3,
                  tol = 1e-6,
                  dampening = c(1,1)) {
    
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
    
    
    ##-----------------------
    #' # initial values 
    ##-----------------------
    fit_start <- gam(list(mu_formula, sigma_formula), data = data, family = gaulss(), method = "ML")
    mu_coefficients <- fit_start$coefficients[1:ncol(mu_X)]^dampening[1] 
    sigma_coefficients <- fit_start$coefficients[-(1:ncol(mu_X))]^dampening[2]
    old_gamma_divergence <- -Inf
    
    
    ##-----------------------
    #' # Run MM-algorithm
    ##-----------------------
    for(r in 1:max_iteration) {
        # weight 
        mu_coefficients_old <- mu_coefficients
        sigma_coefficients_old <- sigma_coefficients
        mu <- as.vector(mu_X %*% mu_coefficients)
        sig <- exp(as.vector(sig_X %*% sigma_coefficients))
        dens <- dnorm(Y, mean = mu, sd = sig)
        ww <- dens^gamma_tuning / mean(dens^gamma_tuning)
        
        # update mu_coefficients
        mu_coefficients <- c( solve(crossprod(mu_X, ww * mu_X / sig^2)) %*% crossprod(mu_X, ww * Y / sig^2) )
        
        # update sigma_coefficients -- Using Fisher scoring here instead of straight up Newton-Raphson
        resid <- Y - as.vector(mu_X %*% mu_coefficients) 
        #Hess <- (-2)*crossprod(X, ww * resid^2 * X / sig^2)
        expFisherinfo <- 2 * crossprod(sig_X, ww * sig_X) 
        Grad <- apply(ww * resid^2 * sig_X / sig^2 - sig_X / (1 + gamma_tuning), 2, sum) 
        #sigma_coefficients <- sigma_coefficients - as.vector(solve(Hess) %*% Grad)
        sigma_coefficients <- sigma_coefficients + as.vector(solve(expFisherinfo) %*% Grad)
        
        # convergence check
        new_gamma_div <- n / gamma_tuning * log(mean(dens^gamma_tuning)) + gamma_tuning / (2 * (1 + gamma_tuning)) * sum(log(sig^2))
        dd <- new_gamma_div - old_gamma_divergence
        if(r %% 10 == 0) 
            message("Iteration: ", r, "\t New gamma divergence: ", round(new_gamma_div, 4),
                    "\t Old gamma divergence: ", round(old_gamma_divergence, 4),
                    "\t Difference: ", round(dd, 4))
        
        old_gamma_divergence <- new_gamma_div
        if(dd < tol) 
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
    
    Z <- (Y - mu)^2 / sig^2 - 1
    a <- ww * (0.5 * gamma_tuning * Z^2 - (Z + 1))
    A <- apply(sig_X * ww * Z, 2, sum)
    Hess_sigma_coefficients <- (-1) * gamma_tuning / (n * S^2) * tcrossprod(A) + (2 / S) * crossprod(sig_X, a * sig_X)
    
    
    ##-----------------------
    #' # Final wrap up 
    ##-----------------------
    C <- ( (1 + gamma_tuning)^(-1/2) * (2*pi*sig^2)^(-gamma_tuning/2) )^(gamma_tuning / (1 + gamma_tuning))
    dens <- dnorm(Y, mean = mu, sd = sig)
    term1 <- dens^(gamma_tuning) * 2*(gamma_tuning * (Y - mu)^2 - sig^2) / ( sig^2 * C) * dens
    term2 <- dens^(2*gamma_tuning) * (Y - mu)^2 / ( sig^2 * C^2) 
    Hscore <- mean(term1 + term2)
    
    dens <- dnorm(Y, mean = mu, sd = sig)
    ww <- dens^gamma_tuning / mean(dens^gamma_tuning)
    
    out <- list(call = match.call(),
                mu_coefficients = as.vector(mu_coefficients), 
                sigma_coefficients = as.vector(sigma_coefficients),
                mu_fitted = mu, 
                sigma_fitted = sig,
                raw_residuals = Y - mu,
                pearson_residuals = (Y - mu) / sig,
                final_weights = ww, 
                gamma_tuning = gamma_tuning,
                gamma_divergence = new_gamma_div, 
                Hscore = Hscore,
                V_mu_coefficients = solve(-Hess_mu_coefficients), 
                V_sigma_coefficients = solve(-Hess_sigma_coefficients))
    names(out$mu_coefficients) <- colnames(mu_X)
    names(out$sigma_coefficients) <- colnames(sig_X)
    rownames(out$V_mu_coefficients) <- colnames(out$V_mu_coefficients) <- colnames(mu_X)
    rownames(out$V_sigma_coefficients) <- colnames(out$V_sigma_coefficients) <- colnames(sig_X)
    
    class(out) <- "gddlm"
    return(out) 
    }
