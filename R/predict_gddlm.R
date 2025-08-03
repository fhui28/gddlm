#' @title Predictions from a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted double linear model produced by [gddlm()] and produces predictions given a new set of values for the model covariates or the original values used for the model fit. Predictions can be accompanied by standard errors and uncertainty intervals based on a large sample normality assumption for the distribution of the model coefficients. The routine can optionally return the matrix by which the model coefficients must be pre-multiplied in order to yield the values of the linear predictor.
#'
#' @param object A fitted `gddlm` object.
#' @param newdata A data frame containing the values of the covariates at which predictions are required. If `newdata` is omitted, predictions are made at the original data used to fit the model.
#' @param predictor Whether component to predict. Options are `mu` or `sigma`.
#' @param type The type of prediction required. The default is on the scale of the linear predictor i.e., the mean if `predictor = "mu"` and the log of the standard deviation if `predictor = "sigma"`. The alternative `response` is on the scale of the component itself i.e., the mean if `predictor = "mu"` and the standard deviation if `predictor = "sigma"`. Finally, the alternative `lpmatrix` returns the matrix by which the model coefficients must be pre-multiplied in order to yield the values of the linear predictor (in this case `se.fit` is ignored).
#' @param se.fit Logical. Should standard errors be returned? Default is \code{FALSE}.
#' @param coverage Numeric. The coverage of the uncertainty intervals. Default is 0.95.
#' @param na.action A function which indicates what should happen when `newdata` contain `NA` values. The default is `na.pass`; see [mgcv::predict.gam()] for more details.].
#' @param ... Currently not used.
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
#' @importFrom mgcv gam gaulss predict.gam
#' @importFrom stats qnorm na.pass
#' @md

predict.gddlm <- function(object,
                          newdata,
                          predictor = "mu",
                          type = "response",
                          se.fit = FALSE,
                          coverage = 0.95,
                          na.action = na.pass,
                          ...) {

    if(!inherits(object, "gddlm"))
        stop("object must be of class \"gddlm\".")

    predictor <- match.arg(predictor, choices = c("mu", "sigma"))
    type <- match.arg(type, choices = c("link", "response", "lpmatrix"))


    #' ## Posthoc fit a double linear model using mgcv solely to get the model matrix
    posthoc_fit_mgcv <- mgcv::gam(list(object$mu_formula,
                                       object$sigma_formula),
                                  data = object$model,
                                  family = gaulss(),
                                  method = "ML")

    get_MM <- mgcv::predict.gam(posthoc_fit_mgcv,
                                newdata = newdata,
                                type = "lpmatrix",
                                na.action = na.action)
    rm(posthoc_fit_mgcv)


    if(predictor == "mu") {
        get_MM <- get_MM[, 1:length(object$mu_coefficients), drop = FALSE]

        if(type == "lpmatrix")
            return(get_MM)

        linpred_prediction <- as.vector(get_MM %*% object$mu_coefficients)
        names(linpred_prediction) <- rownames(get_MM)
        if(type %in% c("link", "response"))
            out <- linpred_prediction

        if(!se.fit)
            return(out)
        if(se.fit) {
            getCov <- object$V_mu_coefficients
            getCov <- 0.5*(getCov + t(getCov))
            v <- rowSums((get_MM %*% getCov) * get_MM)
            names(v) <- rownames(get_MM)
            out <- data.frame(fit = linpred_prediction, se.fit = sqrt(v))

            crit <- qnorm(0.975)

            out$lower <- out$fit - crit * out$se.fit
            out$upper <- out$fit + crit * out$se.fit
            if(type %in% c("link", "response"))
                out <- out

            return(out)
            }
        }

    if(predictor == "sigma") {
        get_MM <- get_MM[, length(object$mu_coefficients) + (1:length(object$sigma_coefficients)), drop = FALSE]

        if(type == "lpmatrix")
            return(get_MM)

        linpred_prediction <- as.vector(get_MM %*% object$sigma_coefficients)
        names(linpred_prediction) <- rownames(get_MM)
        if(type %in% c("link"))
            out <- linpred_prediction
        if(type %in% c("response"))
            out <- exp(linpred_prediction)

        if(!se.fit)
            return(out)
        if(se.fit) {
            getCov <- object$V_sigma_coefficients
            getCov <- 0.5*(getCov + t(getCov))
            v <- rowSums((get_MM %*% getCov) * get_MM)
            names(v) <- rownames(get_MM)
            out <- data.frame(fit = linpred_prediction, se.fit = sqrt(v))

            crit <- qnorm(0.975)

            out$lower <- out$fit - crit * out$se.fit
            out$upper <- out$fit + crit * out$se.fit
            if(type %in% c("link"))
                out <- out
            if(type %in% c("response"))
                out <- exp(out)

            return(out)
            }
        }

    }




