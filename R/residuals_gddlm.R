#' @title Residuals from a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Computes residuals from a fitted `gddlm` object. The residuals can be of different types: `response` corresponds to the raw residuals, `deviance` corresponds to the deviance residuals (note these are slightly ad-hoc as the concept of deviance itself is not clearly defined in the case of gamma divergence; here the deviance residual is computed assuming normally distributed errors, and so **it is equivalent to Pearson residuals**), `pearson` corresponds to the Pearson residuals, and `working` corresponds to the working residuals (which in this case is equivalent to the raw response residuals).
#'
#' @param object A fitted `gddlm` object.
#' @param type The type of residuals to compute. Options are `response`, `deviance`, `pearson`, or `working`.
#' @param ... Currently not used.
#'
#' @return A vector of residuals.
#'
#' @author
#' Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#' @examples
#' \donttest{
#' # See the main `gddlm` help file for examples.
#' }
#'
#' @export
#' @md
residuals.gddlm <- function(object,
                            type = "deviance",
                            ...) {

    if(!inherits(object, "gddlm"))
        stop("model must be of class \"gddlm\".")

    type <- match.arg(type, choices = c("response", "deviance", "pearson", "working"))

    if(type %in% c("response", "working")) {
        res <- object$raw_residuals
        }
    if(type %in% c("pearson", "deviance")) {
        res <- object$pearson_residuals
        }

    return(res)
    }
