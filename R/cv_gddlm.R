#' @title Cross-validation for a gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Performs K-fold cross-validation for gddlm to select the optimal tuning parameter \eqn{gamma} for the gamma divergence. The closer the (chosen) value is to zero, the more the gamma divergence behaves like standard maximum likelihood estimation e.g., the same as what is done in [mgcv::gam()] with the [mgcv::gaulss()] family, or in the [lmls::lmls()] function.
#'
#' @param mu_formula A double-sided formula describing how the response is regressed against the covariates in the mean part of the double linear model.
#' @param sigma_formula A single-sided formula describing how the log of the standard deviation is regressed against the covariates in the double linear model.
#' @param data A data frame containing the response and covariates used in `mu_formula` and `sigma_formula`.
#' @param gamma_seq A numeric vector of candidate values for the tuning parameter \eqn{gamma}.
#' @param nfolds The number of folds to use in the cross-validation. Default is 10.
#' @param fold_seed An optional value to set the random seed for the fold assignment. Default is \code{NULL} i.e., random seed.
#' @param foldid An optional vector of values between 1 and `nfolds` identifying which fold each observation is in. If supplied, then this overrides the `nfolds` arugment. Otherwise it defaults to `NULL`, in which case the folds are randomly assigned.
#' @param ratio_threshold A numeric value that should be greater than 1. To classify how many observations in the test set are to be treated as outliers, the procedure finds the largest ratio of the sorted mean gamma divergence density weights that exceeds this threshold. All observations prior to this are then deemed as outliers in that test fold. Default is 2.
#' @param plotit Should a plot be produced showing the cross-validation error against the candidate values of \eqn{gamma}? Defaults to `FALSE`.
#' @param doparallel Should parallel computing be used to perform the cross-validation? Defaults to `FALSE`.
#' @param num_cores If `do_parallel = TRUE`, then this argument controls the number of cores used. Defaults to \code{NULL}, in which case it is set to `parallel::detectCores() - 2`.
#' @param ... Other arguments to supply to the main [gddlm()] function e.g., `max_iteration` `tolerance`, etc.
#'
#' @details
#' The deviance is used as the measure of predictive accuracy in the cross-validation procedure. The function (can) produces a plot and a return a value for the optimal \eqn{gamma}. Please see manuscript for more details (currently in preparation) as to exactly how the cross-validation procedure works.
#'
#' @return An object of class `cv_gddlm` with the following components:
#' \itemize{
#' \item \code{call} The matched call.
#' \item \code{gamma_seq} The candidate values of \eqn{gamma} supplied by the user.
#' \item \code{gamma_optimal} The optimal value of \eqn{gamma} selected via cross-validation (the value with the lowest median cross-validation error i.e., predictive deviance, across the folds).
#' \item \code{all_CV} A matrix of the cross-validation errors for each fold (rows) and each candidate value of \eqn{gamma} (columns).
#' \item \code{median_CV} The median cross-validation error i.e., predictive deviance. across the folds for each candidate value of \eqn{gamma}.
#' \item \code{foldid} The fold assignment for each observation.
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
#' @importFrom stats median
#' @md
cv_gddlm <- function(mu_formula,
                     sigma_formula,
                     data,
                     gamma_seq = c(0.001, seq(0.05, 0.5, by = 0.05)),
                     nfolds = 10,
                     fold_seed = NULL,
                     foldid = NULL,
                     ratio_threshold = 2,
                     plotit = FALSE,
                     doparallel = FALSE,
                     num_cores = NULL,
                     ...) {

    #' # Checks and balances
    set.seed(fold_seed)

    if(doparallel) {
        if(is.null(num_cores))
            registerDoParallel(cores = detectCores() - 2)
        if(!is.null(num_cores)) {
            registerDoParallel(cores = num_cores)
            }
        }

    if(!is.null(foldid)) {
        message("foldid supplied, overriding nfolds argument. The number of folds is ", length(unique(foldid)), ".")
        fold_index <- foldid
        nfolds <- length(unique(foldid))
        }
    if(is.null(foldid))
        fold_index <- sample(1:nfolds, nrow(data), replace = TRUE)

    #' # Perform cross-validation (non-parallel version)
    if(!doparallel) {
        CV_error <- matrix(NA, nrow = nfolds, ncol = length(gamma_seq))

        for(fold in 1:nfolds) {
            message(paste0("Onto fold ", fold, "..."))
            sel_train <- which(fold_index != fold)
            sel_test <- which(fold_index == fold)

            all_deviance <- all_weights <- matrix(NA, nrow = length(sel_test), ncol = length(gamma_seq))
            for(l in 1:length(gamma_seq)) {
                cw_gamma <- gamma_seq[l]
                fit <- gddlm(mu_formula = mu_formula,
                             sigma_formula = sigma_formula,
                             data = data[sel_train, ],
                             gamma_tuning = cw_gamma,
                             trace = FALSE,
                             ...)

                test_mu <- predict.gddlm(fit,
                                   newdata = data[sel_test,],
                                   predictor = "mu",
                                   type = "response",
                                   se.fit = FALSE)
                test_sigma <- predict.gddlm(fit,
                                      newdata = data[sel_test,],
                                      predictor = "sigma",
                                      type = "response",
                                      se.fit = FALSE)

                fake_gam_mu <- mgcv::gam(mu_formula, data = data[sel_test,], family = gaussian())
                Y <- model.frame(fake_gam_mu)
                Y <- as.vector(model.response(Y))
                rm(fake_gam_mu)

                dens <- dnorm(Y, mean = test_mu, sd = test_sigma)
                ww <- length(sel_test) * dens^cw_gamma / sum(dens^cw_gamma)
                all_weights[, l] <- ww
                all_deviance[,l] <- -2*dnorm(Y, mean = test_mu, sd = test_sigma, log = TRUE)
                rm(Y, test_mu, test_sigma, dens, ww)
                }

            #' ## Determine the mean gamma divergence density weight for each observation in the test test, averaged across all candidate gamma values
            #' Then sort these weights, and compute the ratio of the current average weight on the previous average weight
            sel_non_outliers <- sort(rowMeans(all_weights))[-1]/sort(rowMeans(all_weights))[-nrow(all_weights)]

            #' ## If none of the ratios are greater than a pre-specified threshold, then use all observations in the test set to compute the CV error for this fold
            #' Otherwise, find the largest ratio exceeding the pre-specified threshold, and declare all observations prior to this as being outliers. The remaining observations are then used to compute the CV error for this fold.
            #' The intuition is that if there is a large jump in the average weights, then this indicates that the observations prior to this jump are outliers (because they weights are relative small compared to the remaining observations).
            if(all(sel_non_outliers < ratio_threshold))
                final_selection <- 1:nrow(all_weights)
            if(any(sel_non_outliers > ratio_threshold)) {
                sel_non_outliers <- which.max(sel_non_outliers)
                final_selection <- order(rowMeans(all_weights))[-(1:sel_non_outliers)]
                }

            #' CV error is the mean deviance for the selected observations in the test set
            final_deviance <- colMeans(all_deviance[final_selection,])
            CV_error[fold, ] <- final_deviance
            rm(all_weights, all_deviance, sel_non_outliers, final_selection, final_deviance, fit, cw_gamma)
            }
        }

        #' # Perform cross-validation (parallel version)
    if(doparallel) {
        compute_cv_error_fn <- function(fold) {
            sel_train <- which(fold_index != fold)
            sel_test <- which(fold_index == fold)

            all_deviance <- all_weights <- matrix(NA, nrow = length(sel_test), ncol = length(gamma_seq))
            for(l in 1:length(gamma_seq)) {
                cw_gamma <- gamma_seq[l]
                fit <- gddlm(mu_formula = mu_formula,
                             sigma_formula = sigma_formula,
                             data = data[sel_train, ],
                             gamma_tuning = cw_gamma,
                             trace = FALSE,
                             ...)

                test_mu <- predict.gddlm(fit,
                                   newdata = data[sel_test,],
                                   predictor = "mu",
                                   type = "response",
                                   se.fit = FALSE)
                test_sigma <- predict.gddlm(fit,
                                      newdata = data[sel_test,],
                                      predictor = "sigma",
                                      type = "response",
                                      se.fit = FALSE)

                fake_gam_mu <- mgcv::gam(mu_formula, data = data[sel_test,], family = gaussian())
                Y <- model.frame(fake_gam_mu)
                Y <- as.vector(model.response(Y))
                rm(fake_gam_mu)

                dens <- dnorm(Y, mean = test_mu, sd = test_sigma)
                ww <- length(sel_test) * dens^cw_gamma / sum(dens^cw_gamma)
                all_weights[, l] <- ww
                all_deviance[,l] <- -2*dnorm(Y, mean = test_mu, sd = test_sigma, log = TRUE)
                rm(Y, test_mu, test_sigma, dens, ww)
                }

            sel_non_outliers <- sort(rowMeans(all_weights))[-1]/sort(rowMeans(all_weights))[-nrow(all_weights)]
            if(all(sel_non_outliers < ratio_threshold))
                final_selection <- 1:nrow(all_weights)
            if(any(sel_non_outliers > ratio_threshold)) {
                sel_non_outliers <- which.max(sel_non_outliers)
                final_selection <- order(rowMeans(all_weights))[-(1:sel_non_outliers)]
                }

            #' CV error is the mean deviance for the selected observations in the test set
            return(colMeans(all_deviance[final_selection,]))
            }

        CV_error <- foreach(fold = 1:nfolds, .combine = rbind) %dopar% {
            compute_cv_error_fn(fold)
            }
        }

    set.seed(NULL)
    out <- list(call = match.call(),
                gamma_seq = gamma_seq,
                gamma_optimal = gamma_seq[which.min(apply(CV_error, 2, median))],
                all_CV = CV_error,
                median_CV = apply(CV_error, 2, median),
                foldid = fold_index)
    if(plotit) {
        p <- ggplot(data.frame(gamma_seq = out$gamma_seq,
                               median_cv = out$median_CV),
                    aes(x = .data$gamma_seq, y = .data$median_cv)) +
            geom_point() +
            geom_line() +
            labs(x = expression(gamma), y = "Median CV Error") +
            theme_bw()
        print(p)
        }

    class(out) <- "cv_gddlm"
    return(out)
    }




