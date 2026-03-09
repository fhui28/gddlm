#' @title Model diagnostics plots for gddlm models
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is basically a wrapper around the [gratia::appraise()] function, but designed to work with `gddlm` models. The main difference is that observations (residuals, quantiles etc...) in the diagnostic plots are colored by the final weights from the `gddlm` model.
#'
#' @param model A fitted `gddlm` model.
#' @param use_worm Logical. Should the worm plot be drawn in place of the QQ plot? Default is `FALSE`.
#' @param n_uniform Integer. The number of times to randomize uniform quantiles in the direct computation method or QQ plots.
#' @param n_simulate Integer. The number of times to simulate from the estimated model in the simulation-based method for QQ plots.
#' @param method Character. The method to use for to generate theoretical quantiles. Options are "uniform" or "simulate".
#' @param type Character. The type of residuals to use. Options are "deviance" or "pearson" residuals.
#' @param n_bins Character. The number of bins to use in the histogram of residuals. See [gratia::appraise()] for more details.
#' @param ncol Integer. The number of columns in the grid of plots. If `NULL`, the number of columns is determined automatically.
#' @param nrow Integer. The number of rows in the grid of plots. If `NULL`, the number of rows is determined automatically.
#' @param guides Character. The guide to use for the color scale. See [gratia::appraise()] for more details
#' @param level Numeric. The confidence level for the reference intervals in the QQ plot. Only used if `method = "simulate"`.
#' @param ci_col Character. The color of the reference interval lines in the QQ plot. Only used if `method = "simulate"`.
#' @param ci_alpha Numeric. The alpha transparency level of the reference interval lines in the QQ plot. Only used if `method = "simulate"`.
#' @param line_col Character. The color of the reference 1:1 line in the QQ plot, and the reference line in the residuals vs. linear predictor plot.
#' @param ... Currently not used.
#'
#'
#' @note
#' Acknowledgments to Gavin Simpson and the `gratia` package for the original [gratia::appraise()] function.
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
#' @importFrom gratia appraise
#' @importFrom colorspace scale_color_continuous_sequential
#' @importFrom ggplot2 geom_point theme_bw aes_string
#' @importFrom mgcv gam gaulss
#' @md

appraise.gddlm <- function(model,
                           use_worm = FALSE,
                           n_uniform = 10,
                           n_simulate = 50,
                           method = "uniform",
                           type = "deviance",
                           n_bins = "sturges",
                           ncol = NULL,
                           nrow = NULL,
                           guides = "keep",
                           level = 0.9,
                           ci_col = "black",
                           ci_alpha = 0.2,
                           line_col = "red",
                           ...) {

    if(!inherits(model, "gddlm"))
        stop("model must be of class \"gddlm\".")

    ##-----------------------
    #' ## Duplicate GAM model but remove as many things as you do not need (mainly done so that gratia does not take in the wrong information)
    ##-----------------------
    posthoc_mgcv_fit_fn <- function(model) {
        dofit <- mgcv::gam(list(model$mu_formula,
                                model$sigma_formula),
                           data = model$model,
                           family = gaulss(),
                           method = "ML")
        return(dofit)
    }
    posthoc_mgcv_fit <- posthoc_mgcv_fit_fn(model = model)

    posthoc_mgcv_fit$aic <- NULL
    posthoc_mgcv_fit$assign <- NULL
    posthoc_mgcv_fit$boundary <- NULL
    posthoc_mgcv_fit$call <- NULL
    posthoc_mgcv_fit$control <- NULL
    posthoc_mgcv_fit$converged <- NULL
    posthoc_mgcv_fit$db.drho <- NULL
    posthoc_mgcv_fit$deviance <- NULL
    posthoc_mgcv_fit$df.null <- NULL
    posthoc_mgcv_fit$df.residual <- NULL
    posthoc_mgcv_fit$edf <- NULL
    posthoc_mgcv_fit$edf1 <- NULL
    posthoc_mgcv_fit$edf2 <- NULL
    posthoc_mgcv_fit$full.sp <- NULL
    posthoc_mgcv_fit$F <- NULL
    posthoc_mgcv_fit$gcv.ubre <- NULL
    posthoc_mgcv_fit$hat <- NULL
    posthoc_mgcv_fit$iter <- NULL
    posthoc_mgcv_fit$method <- NULL
    posthoc_mgcv_fit$mgcv.conv <- NULL
    posthoc_mgcv_fit$min.edf <- NULL
    posthoc_mgcv_fit$nsdf <- NULL
    posthoc_mgcv_fit$null.deviance <- NULL
    posthoc_mgcv_fit$optimizer <- NULL
    posthoc_mgcv_fit$outer.info <- NULL
    posthoc_mgcv_fit$paraPen <- NULL
    posthoc_mgcv_fit$pterms <- NULL
    posthoc_mgcv_fit$R <- NULL
    posthoc_mgcv_fit$rank <- NULL
    posthoc_mgcv_fit$reml.scale <- NULL
    posthoc_mgcv_fit$residuals <- NULL
    posthoc_mgcv_fit$rV <- NULL
    posthoc_mgcv_fit$smooth <- NULL
    posthoc_mgcv_fit$var.summary <- NULL
    posthoc_mgcv_fit$Ve <- NULL
    posthoc_mgcv_fit$Vp <- NULL
    posthoc_mgcv_fit$Vc <- NULL
    posthoc_mgcv_fit$weights <- NULL


    ##-----------------------
    #' ## Modify relevant items and pass to gratia::appraise
    ##-----------------------
    posthoc_mgcv_fit$coefficients <- c(model$mu_coefficients, model$sigma_coefficients)
    names(posthoc_mgcv_fit$coefficients)[1:length(model$mu_coefficients)] <- paste0("mu_", names(model$mu_coefficients))
    names(posthoc_mgcv_fit$coefficients)[-(1:length(model$mu_coefficients))] <- paste0("sigma_", names(model$sigma_coefficients))
    posthoc_mgcv_fit$Vp <- as.matrix(Matrix::bdiag(model$V_mu_coefficients, model$V_sigma_coefficients))
    posthoc_mgcv_fit$Ve <- as.matrix(Matrix::bdiag(model$V_mu_coefficients, model$V_sigma_coefficients))
    posthoc_mgcv_fit$Vc <- as.matrix(Matrix::bdiag(model$V_mu_coefficients, model$V_sigma_coefficients))

    posthoc_mgcv_fit$fitted.values <- cbind(model$mu_fitted, 1/model$sigma_fitted) #' See help file for fitted values with gaulss
    posthoc_mgcv_fit$linear.predictors <- cbind(posthoc_mgcv_fit$fitted.values[,1], #' See help file for fitted values with gaulss
                                                log(1/posthoc_mgcv_fit$fitted.values[,2]))
    posthoc_mgcv_fit$residuals <- residuals.gddlm(model, type = type)

    make_plots <- gratia::appraise(posthoc_mgcv_fit,
                                   method = method,
                                   type = "deviance", # Since Pearson and deviance residuals are defined the same way for the postdoc double linear model fit
                                   use_worm = use_worm,
                                   n_uniform = n_uniform,
                                   n_simulate = n_simulate,
                                   n_bins = n_bins,
                                   ncol = ncol,
                                   nrow = nrow,
                                   guides = guides,
                                   level = level,
                                   ci_col = ci_col,
                                   ci_alpha = ci_alpha,
                                   line_col = line_col)

    make_plots[[1]]$data$weights <- model$final_weights[order(posthoc_mgcv_fit$residuals)]
    make_plots[[1]]$layers[[2]] <- NULL
    make_plots[[1]] <- make_plots[[1]] +
        ggplot2::geom_point(aes_string(x = 'theoretical', y = 'residuals', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    make_plots[[2]]$data$weights <- model$final_weights
    make_plots[[2]]$layers[[2]] <- NULL
    make_plots[[2]] <- make_plots[[2]] +
        ggplot2::geom_point(aes_string(x = 'eta', y = 'residuals', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    make_plots[[3]] <- make_plots[[3]] +
        ggplot2::theme_bw()

    make_plots[[4]]$data$weights <- model$final_weights
    make_plots[[4]]$layers[[2]] <- NULL
    make_plots[[4]] <- make_plots[[4]] +
        ggplot2::geom_point(aes_string(x = 'fitted', y = 'observed', color = 'weights')) +
        colorspace::scale_color_continuous_sequential(palette = "Plasma") +
        ggplot2::theme_bw()

    print(make_plots)
}


