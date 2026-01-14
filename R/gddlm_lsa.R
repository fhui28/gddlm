#' @title Fast robust model selection for double linear models using gamma divergence.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' After applying [gddlm()] to fit a robust double linear model via gamma divergence, this function constructs a regularization path for the model coefficients based on applying the [broken adaptive ridge or BAR penalty](https://www3.stat.sinica.edu.tw/statistica/J30N2/J30N223/J30N223.html) to a [least squares (quadratic) approximation](https://www.tandfonline.com/doi/abs/10.1198/016214507000000509) of the gamma divergence objective function. This approach is a relatively computationally efficient approach to performing joint selection on the coefficients characterizing the mean (mu) and the log of the standard deviation (sigma) of the double linear model.
#'
#' @param object An object of class `gddlm`.
#' @param nlambda The number of tuning parameters values.
#' @param lambda_min_ratio The smallest value for the tuning parameter lambda, as a fraction of the maximum internally-derived lambda value i.e., a value such that penalization leads to zeros for all penalized estimates.
#' @param lambda A user-supplied tuning parameter sequence. Note this is usually not supplied, as it is standard to have the function itself compute its own lambda sequence based on `nlambda` and `lambda_min_ratio`.
#' @param min_df The minimum number of non-zero coefficients, excluding the two intercepts in the mu and sigma components, allowed in the double linear model. This is useful to supply if, when the function tries to find a appropriate lambda sequence, the user wants the maximum value of lambda to not necessarily shrink to all the penalized estimates to zero. Defaults to zero, and can take a maximum value given by `sum(prior_weights != 0)`.
#' @param num_cores The regularization path is constructed via parallel computing, and this argument controls the number of cores used. Defaults to `NULL`, in which case it is set to `parallel::detectCores() - 2`.
#' @param prior_weights A vector of prior weights for the BAR penalty, with length equal to the number of parameters in the model (i.e., length of `c(object$mu_coefficients, object$sigma_coefficients)`). Default is `NULL`, in which case all penalized parameters are given equal weight of 1, and the two intercepts (one for mu and one for sigma) are given weight zero (i.e., not penalized). It is strongly recommended to only set elements of `prior_weights` to 0 or 1, as the BAR penalty usually does not need adaptive/fractional weights for selection purposes.
#' @param control A list containing the following elements:
#' \itemize{
#' \item{max_iteration: }{the maximum number of iterations in the BAR penalty algorithm.}
#' \item{tol: }{the convergence criterion for the BAR penalty algorithm; the squared norm of the difference between estimates from successive iterations must be smaller than this value.}
#' \item{round_eps: }{a tolerance to round values to zero; this is technically not needed as the BAR penalty will produce exactly zero estimates up to machine error, but is included anyway.}
#' }
#'
#'
#' @details
#' See associated manuscript (currently in revision).
#'
#' In terms of model selection, the function computes the following information criteria at each point on the regularization path:
#'
#' - Akaike information criterion (AIC) with model complexity penalty of 2;
#' - Bayesian information criterion (BIC) with model complexity penalty of `log(N)` where `N` is the sample size;
#' - Extended Bayesian information criterion (EBIC) with model complexity penalty of `log(N) + log(penalized parameters)`;
#' - Extended Bayesian information criterion (EBIC2) with model complexity penalty of `log(N) + 2 * log(penalized parameters)`.
#'
#'
#' @return An object of class `gddlm_lsa` with the following elements (and not necessarily in the order below):
#' \item{call:}{The function call.}
#' \item{lambda:}{The actual sequence of tuning parameter, lambda, values used.}
#' \item{mu_coefficients_path:}{A sparse matrix showing the estimated parameters corresponding to the mean (mu) component at each point on the regularization path. The number of rows is equal `length(object$mu_coefficients)`, while the number of columns is equal to `length(lambda)`.}
#' \item{sigma_coefficients_path:}{A sparse matrix showing the estimated parameters corresponding to the log standard deviation (sigma) component at each point on the regularization path. The number of rows is equal `length(object$sigma_coefficients)`, while the number of columns is equal to `length(lambda)`.}
#' \item{information:}{A data frame containing various information about the regularization path, including:
#' The sequence of lambda values used;
#' The value of the (unpenalized) least squares approximation statistic at each point on the regularization path;
#' The degrees of freedom (number of non-zero estimated parameters) at each point on the regularization path;
#' Various information criteria at each point on the regularization path (see the details above for more information of these).
#' }
#' \item{time_taken:}{The time taken to construct the regularization path.}
#'
#' @author Francis K.C. Hui <fhui28@gmail.com>
#'
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(GJRM)
#' library(gamlss.lasso)
#' library(mgcv)
#' library(mvtnorm)
#'
#' # Simulation setting and data generation
#' set.seed(012026)
#' n <- 500
#' p <- 25
#' H <- abs(outer(1:p, 1:p, "-"))
#' H <- 0.75^H
#' X <- cbind(1, rmvnorm(n, sigma = H))
#' Beta_true <- c(-1, -1, 1, -0.5, 0.5, 0, 0, 0.25, rep(0, p - 7))
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
#'
#' # Construct regularization path via least squares approximation with BAR penalty
#' fit_gddlm_choosegamma_lsa <- gddlm_lsa(object = fit_gddglm_choosegamma)
#'
#' ## Example of extracting estimates based minimizing BIC
#' estimates(fit_gddlm_choosegamma_lsa, lambda = "BIC")
#'
#'
#' # Fitting some other double linear models and attempt variable selection using this...
#' ## Non-robust double linear model with selection using the adaptive lasso
#' dat <- data.frame(y = Y, X[,-1])
#' fit_gamlsslasso <- gamlss::gamlss(y ~ gnet(x.vars = names(dat)[-1]),
#' sigma.fo = ~ gnet(x.vars = names(dat)[-1]),
#' data = dat,
#' family = NO,
#' i.control = glim.control(cyc=1, bf.cyc=1))
#'
#'
#' ## Psi-divergence. Note this takes the longest by far...
#' ## Then perform varable based on simple Wald-tests from the saturated model
#' fit_psidiv <- GJRM::gamlss(list(make_mu_formula, make_sigma_formula),
#' data = data.frame(response = Y, X[,-1]),
#' family = "N",
#' robust = TRUE)
#'
#'
#' ## Assess performance
#' psidiv_coef <- fit_psidiv$coefficients[1:(p+1)]
#' psidiv_coef[abs(fit_psidiv$coefficients[1:(p+1)] / sqrt(diag(fit_psidiv$Vb))[1:(p+1)]) < qnorm(0.975)] <- 0
#' Beta_est = data.frame(
#' true = Beta_true,
#' gddlm = estimates(fit_gddlm_choosegamma_lsa, lambda = "BIC")$mu_coefficients,
#' psidiv = psidiv_coef,
#' gamlss_lasso = c(NA, tail(getSmo(fit_gamlsslasso, "mu") ,1)[[1]]$beta)
#' )
#' rm(psidiv_coef)
#' Beta_est
#'
#'
#' psidiv_coef <- fit_psidiv$coefficients[-(1:(p+1))]
#' psidiv_coef[abs(fit_psidiv$coefficients[-(1:(p+1))] / sqrt(diag(fit_psidiv$Vb))[-(1:(p+1))])] <- 0
#' Delta_est = data.frame(
#' true = Delta_true,
#' gddlm = estimates(fit_gddlm_choosegamma_lsa, lambda = "BIC")$sigma_coefficients,
#' psidiv = psidiv_coef,
#' gamlss_lasso = c(NA, tail(getSmo(fit_gamlsslasso, "sigma") ,1)[[1]]$beta)
#' )
#' rm(psidiv_coef)
#' Delta_est
#'
#' }
#'
#'
#' @export gddlm_lsa
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom parallel detectCores
#' @import Matrix
#' @md

gddlm_lsa <- function(object,
                      nlambda = 100,
                      lambda_min_ratio = 1e-6,
                      lambda = NULL,
                      min_df = 0,
                      num_cores = NULL,
                      prior_weights = NULL,
                      control = list(max_iteration = 1000, tol = 1e-6, round_eps = 1e-6)) {

    ##-------------------
    #' # Checks and balances
    ##-------------------
    if(!inherits(object, "gddlm"))
        stop("object must be of class 'gddlm'.")
    if(is.null(num_cores))
        num_cores <- detectCores()-2
    registerDoParallel(num_cores)
    message("Using ", num_cores, " cores in construct regularization path.")

    if(is.null(control$max_iteration))
        control$max_iteration <- 100
    if(is.null(control$tol))
        control$tol <- 1e-6
    if(is.null(control$round_eps))
        control$round_eps <- 1e-6

    psihat <- c(object$mu_coefficients, object$sigma_coefficients)
    coefficient_names <- names(psihat)
    if(is.null(prior_weights)) {
        prior_weights <- rep(1, length(psihat))
        prior_weights[c(1, length(object$mu_coefficients) + 1)] <- 0  # Do not penalize intercepts
        }
    if(length(prior_weights) != length(psihat))
        stop("Length of prior_weights must be equal to the number of parameters in the model.")
    if(min_df > sum(prior_weights != 0))
        stop("min_df cannot be larger than the number of penalized parameters in the model.")


    #' ## Construct bread/Hessian matrix
    bread_matrix <- bdiag(solve(object$V_mu_coefficients),
                          solve(object$V_sigma_coefficients))


    ##-----------------------
    ## Attempt to determine a minimum and maximum lambda, if required, for BAR penalty
    ##-----------------------
    if(!is.null(lambda))
        lambdaseq <- lambda
    if(is.null(lambda)) {
        message("Finding an appropriate sequence of lambda values...")
        lambda_max <- 0.05

        fit_togetmaxlambda <- .BAR_singlelambda(psihat = psihat,
                                                bread_matrix = bread_matrix,
                                                prior_weights = prior_weights,
                                                lambda = lambda_max,
                                                control = control)

        while(sum(fit_togetmaxlambda[prior_weights != 0] != 0) > min_df) {
            lambda_max <- lambda_max*1.1
            fit_togetmaxlambda <- .BAR_singlelambda(psihat = psihat,
                                                    bread_matrix = bread_matrix,
                                                    prior_weights = prior_weights,
                                                    lambda = lambda_max,
                                                    control = control)
            }

        lambdaseq <- lseq(lambda_max,
                          lambda_max*lambda_min_ratio,
                          length = nlambda,
                          decreasing = TRUE)
        rm(fit_togetmaxlambda)
        }

    ##-----------------------
    #' ## Construct regularization path for BAR penalty
    ##-----------------------
    tic <- proc.time()
    message("Constructing regularization path...")
    cwfit_fn <- function(l) {
        cwfit <- .BAR_singlelambda(psihat = psihat,
                                   bread_matrix = bread_matrix,
                                   prior_weights = prior_weights,
                                   lambda = lambdaseq[l],
                                   control = control)

        out <- list(
            cw_coefficients = cwfit,
            df_path = sum(cwfit != 0),
            LSA_statistic = as.vector(t(as.vector(cwfit - psihat)) %*% bread_matrix %*% as.vector(cwfit - psihat))
            )

        out$AIC <- out$LSA_statistic + 2*out$df_path
        out$BIC <- out$LSA_statistic + log(nrow(object$model))*out$df_path
        out$EBIC <- out$LSA_statistic + (log(nrow(object$model)) + log(sum(prior_weights != 0)))*out$df_path #' Note log(P) counts only the number of penalized parameters
        out$EBIC2 <- out$LSA_statistic + (log(nrow(object$model)) + 2 * log(sum(prior_weights != 0)))*out$df_path
        return(out)
        }

    allfits <- foreach(l = 1:length(lambdaseq)) %dopar% cwfit_fn(l = l)
    toc <- proc.time()

    coefficients_path <- Matrix::Matrix(sapply(allfits, function(x) x$cw_coefficients), sparse = TRUE)
    rownames(coefficients_path) <- names(psihat)
    df_path <- sapply(allfits, function(x) x$df_path)
    LSQ_statistics <- sapply(allfits, function(x) x$LSA_statistic)
    AIC <- sapply(allfits, function(x) x$AIC)
    BIC <- sapply(allfits, function(x) x$BIC)
    EBIC <- sapply(allfits, function(x) x$EBIC)
    EBIC2 <- sapply(allfits, function(x) x$EBIC2)

    out <- list(lambda = lambdaseq,
                mu_coefficients_path = coefficients_path[1:length(object$mu_coefficients), ],
                sigma_coefficients_path = coefficients_path[length(object$mu_coefficients) + 1:length(object$sigma_coefficients), ],
                information = data.frame(lambda = lambdaseq,
                                         LSQ_statistic = LSQ_statistics,
                                         df = df_path,
                                         AIC,
                                         BIC,
                                         EBIC,
                                         EBIC2),
                time_taken = toc - tic)

    ##-----------------------
    ## Final touch up
    ##-----------------------
    rownames(out$mu_coefficients_path) <- names(object$mu_coefficients)
    rownames(out$sigma_coefficients_path) <- names(object$sigma_coefficients)

    out$call <- match.call()
    class(out) <- "gddlm_lsa"
    return(out)
    }



#' @title Hidden internal function for applying the BAR penalty applied to a least squares approximation of a GDDLM
#' @description Adapted from hidden internal function in (https://github.com/fhui28/COQUE)
#' @return A vector of penalized estimates with the same length as \code{psihat}, with potentially one or more elements shrunk to zero.
#' @importFrom methods as
#' @noRd
#' @noMd

.BAR_singlelambda <- function(psihat,
                              bread_matrix,
                              prior_weights,
                              lambda,
                              control) {

    ##-------------------------
    #' ## Attempting the analogue of equation 4 in [https://doi.org/10.1016/j.jmva.2018.08.007]
    ##-------------------------
    Dbar <- Diagonal(x = prior_weights * lambda, n = length(psihat))
    err <- Inf
    counter <- 0

    cw_psi <- psihat
    while(counter < control$max_iteration & err > control$tol) {
        GammaMatrix_psi <- diag(x = cw_psi, nrow = length(psihat))
        diag(GammaMatrix_psi)[prior_weights == 0] <- 1
        GammaMatrix_psi <- as(GammaMatrix_psi, "sparseMatrix")

        new_psi <- GammaMatrix_psi %*% solve(GammaMatrix_psi %*% bread_matrix %*% GammaMatrix_psi + Dbar) %*% GammaMatrix_psi %*% bread_matrix %*% psihat
        new_psi <- as.vector(new_psi)
        names(new_psi) <- names(psihat)

        err <- sum((new_psi - cw_psi)^2)
        cw_psi <- new_psi
        counter <- counter + 1
        }

    #' Note the convergence properties of the algorithm does basically ensure exact zeros can be achieved (up to machine error); see [https://doi.org/10.1155/2016/3456153]. However once can also round to zero if interested. In this case, I will round based on eps = 1e-6
    new_psi[abs(new_psi) < control$round_eps] <- 0
    return(new_psi)
    }


