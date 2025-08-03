#' @title Calculate covariance matrix for a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Extracts the covariance (inverse Fisher information) matrix corresponding to either the mean or the log standard deviation part of the fitted robust double linear model.
#'
#'
#' @param object A fitted `gddlm` object.
#' @param type Which component to extract. Options are `mu`, or `sigma`.
#' @param ... Currently not used.
#'
#' @return A vector of fitted values
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
#' @md
vcov.gddlm <- function(object,
                       type = "mu",
                       ...) {

    if(!inherits(object, "gddlm"))
        stop("model must be of class \"gddlm\".")

    type <- match.arg(type, choices = c("mu", "sigma"))

    if(type %in% c("mu")) {
        res <- object$V_mu_coefficients
    }
    if(type %in% c("sigma")) {
        res <- object$V_sigma_coefficients
    }

    return(res)
}


#' @rdname fitted.gddlm
#' @export
fitted.values.gddlm <- fitted.gddlm
