#' @title Choosing gamma for a gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Performs a pre-fit using `gddlm` with a relatively large value of `gamma_tuning` to first identify a small set of non-contaminated (clean) observations. Then, for each candidate value of `gamma_seq`, the function fits a `gddlm` model on the rest of the data and computes the mean deviance on the clean observations. To assess variability around this, bootstrap samples of the clean observations and the mean deviance evaluated on for each bootstrap sample. The optimal value of \eqn{\gamma} is then selected as the one with the lowest median (or mean) cross-validation error across the bootstrap samples.
#'
#' @param mu_formula A double-sided formula describing how the response is regressed against the covariates in the mean part of the double linear model.
#' @param sigma_formula A single-sided formula describing how the log of the standard deviation is regressed against the covariates in the double linear model.
#' @param data A data frame containing the response and covariates used in `mu_formula` and `sigma_formula`.
#' @param gamma_seq A numeric vector of candidate values for the tuning parameter \eqn{\gamma}.
#' @param prefit_gamma The value of \eqn{\gamma} used in the pre-fit to identify the clean observations. Default is 0.5, which from personal experience is sufficient large in many datasets and seems to work well in practice.
#' @param cleanobs_quantile The quantile of the weights from the pre-fit used to identify the clean observations. Default is 0.95, which means that the top 5% of observations with the largest weights are considered as clean observations.
#' @param num_boot The number of bootstrap samples to use. Default is 50.
#' @param plotit Should a plot be produced showing the cross-validation error against the candidate values of \eqn{\gamma}? Defaults to `FALSE`.
#' @param doparallel Should parallel computing be used to perform the cross-validation? Defaults to `FALSE`.
#' @param num_cores If `do_parallel = TRUE`, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to `parallel::detectCores() - 2`.
#' @param ... Other arguments to supply to the main [gddlm()] function e.g., `control` etc.
#'
#' @details
#' The deviance is used as the measure of predictive accuracy in the procedure. The function (can) produces a plot and a return a value for the optimal \eqn{\gamma}.
#'
#' Please see manuscript for more details (currently in preparation) as to exactly how the cross-validation procedure works.
#'
#' @return An object of class `choosegamma_gddlm` with the following components:
#' \itemize{
#' \item \code{call} The matched call.
#' \item \code{gamma_seq} Same as the input argument `gamma_seq`.
#' \item \code{gamma_optimal_mean} The optimal value of \eqn{\gamma} selected using the bootstrap, based on the with the lowest mean predictive deviance, across the bootstrap samples.
#' \item \code{gamma_optimal_median} The optimal value of \eqn{\gamma} selected using the bootstrap, based on the with the lowest median predictive deviance, across the bootstrap samples.
#' \item \code{all_bootstrap_deviance} A matrix of the bootstrap predictive deviances for each bootstrap dataset (rows) and each candidate value of \eqn{\gamma} (columns).
#' }
#'
#' @author
#' Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#' @examples
#' \dontrun{
#' See the main `gddlm` help file for examples.
#' }
#'
#' @export
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom ggplot2 ggplot geom_point geom_line theme_bw labs aes
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom rlang .data
#' @importFrom stats median quantile
#' @md
choosegamma_gddlm <- function(mu_formula,
                              sigma_formula,
                              data,
                              gamma_seq = c(0.001, seq(0.05, 0.5, by = 0.05)),
                              prefit_gamma = 0.5,
                              cleanobs_quantile = 0.95,
                              num_boot = 50,
                              plotit = FALSE,
                              doparallel = FALSE,
                              num_cores = NULL,
                              ...) {

    ##-----------------------
    #' # Checks and balances
    ##-----------------------
    if(doparallel) {
        if(is.null(num_cores))
            num_cores <- detectCores() - 2
        registerDoParallel(cores = num_cores)
        message("Using ", num_cores, " cores in construct regularization path.")
        }

    if(cleanobs_quantile < 0 | cleanobs_quantile > 1)
        stop("cleanobs_quantile must be between 0 and 1.")


    ##-----------------------
    #' # Perform GDDLM prefit and identify non-contaminated observations
    ##-----------------------
    prefit <- gddlm(mu_formula = mu_formula,
                 sigma_formula = sigma_formula,
                 data = data,
                 gamma_tuning = prefit_gamma,
                 ...)
    sel_cleanobs <- which(prefit$final_weights > quantile(prefit$final_weights, cleanobs_quantile))
    sel_cleanobs_boot <- replicate(num_boot, sample(sel_cleanobs, size = length(sel_cleanobs), replace = TRUE))
    fake_gam_mu <- mgcv::gam(mu_formula, data = data, family = gaussian())
    Y <- model.frame(fake_gam_mu)
    Y <- as.vector(model.response(Y))
    rm(fake_gam_mu)


    ##-----------------------
    #' # Perform fitting on rest of the data and compute the (bootstrap) error on the clean dataset for each candidate value of gamma
    ##-----------------------
    if(!doparallel) {
        all_boot_deviance <- matrix(NA, nrow = num_boot, ncol = length(gamma_seq))

        for(l in 1:length(gamma_seq)) {
            cw_gamma <- gamma_seq[l]
            fit_rest <- gddlm(mu_formula = mu_formula,
                              sigma_formula = sigma_formula,
                              data = data[-sel_cleanobs, ],
                              gamma_tuning = cw_gamma,
                              ...)

            test_mu <- predict.gddlm(fit_rest,
                                     newdata = data[sel_cleanobs,],
                                     predictor = "mu",
                                     type = "response",
                                     se.fit = FALSE)
            test_mu <- matrix(test_mu[as.character(sel_cleanobs_boot)], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot))
            test_sigma <- predict.gddlm(fit_rest,
                                        newdata = data[sel_cleanobs,],
                                        predictor = "sigma",
                                        type = "response",
                                        se.fit = FALSE)
            test_sigma <- matrix(test_sigma[as.character(sel_cleanobs_boot)], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot))
            all_boot_deviance[,l] <- -2 * apply(dnorm(matrix(Y[sel_cleanobs_boot], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot)),
                                                mean = test_mu,
                                                sd = test_sigma, log = TRUE), 2, mean, na.rm = TRUE)
            rm(test_mu, test_sigma)
            }
        }


    ##-----------------------
    #' # Perform fitting on rest of the data and compute the (bootstrap) error on the clean dataset for each candidate value of gamma (parallel version)
    ##-----------------------
    if(doparallel) {
        compute_boot_deviance_fn <- function(l) {
            cw_gamma <- gamma_seq[l]
            fit_rest <- gddlm(mu_formula = mu_formula,
                              sigma_formula = sigma_formula,
                              data = data[-sel_cleanobs, ],
                              gamma_tuning = cw_gamma,
                              ...)

            test_mu <- predict.gddlm(fit_rest,
                                     newdata = data[sel_cleanobs,],
                                     predictor = "mu",
                                     type = "response",
                                     se.fit = FALSE)
            test_mu <- matrix(test_mu[as.character(sel_cleanobs_boot)], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot))
            test_sigma <- predict.gddlm(fit_rest,
                                        newdata = data[sel_cleanobs,],
                                        predictor = "sigma",
                                        type = "response",
                                        se.fit = FALSE)
            test_sigma <- matrix(test_sigma[as.character(sel_cleanobs_boot)], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot))
            boot_deviance <- -2 * apply(dnorm(matrix(Y[sel_cleanobs_boot], nrow = nrow(sel_cleanobs_boot), ncol = ncol(sel_cleanobs_boot)),
                                             mean = test_mu,
                                             sd = test_sigma, log = TRUE), 2, mean, na.rm = TRUE)
            rm(test_mu, test_sigma)

            return(boot_deviance)
            }

        all_boot_deviance <- foreach(l = 1:length(gamma_seq), .combine = cbind) %dopar% { compute_boot_deviance_fn(l) }
        }


    ##-----------------------
    #' Final output
    ##-----------------------
    set.seed(NULL)
    out <- list(call = match.call(),
                gamma_seq = gamma_seq,
                gamma_optimal_mean = gamma_seq[which.min(apply(all_boot_deviance, 2, mean))],
                gamma_optimal_median = gamma_seq[which.min(apply(all_boot_deviance, 2, median))],
                all_bootstrap_deviance = all_boot_deviance)
    if(plotit) {
        p <- ggplot(data.frame(gamma_seq = out$gamma_seq,
                               boot_deviance_median = apply(all_boot_deviance, 2, median)),
                    aes(x = .data$gamma_seq, y = .data$boot_deviance_median)) +
            geom_point() +
            geom_line() +
            labs(x = expression(gamma), y = "Median CV Error") +
            theme_bw()
        print(p)
        }

    class(out) <- "choosegamma_gddlm"
    return(out)
    }




