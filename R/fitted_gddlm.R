#' @title Extract fitted values from a gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Extracts fitted values corresponding to either the mean or the log standard deviation part of the fitted robust double linear model. 
#'
#'
#' @param object A fitted `gddlm` object.
#' @param type Which component to extract. Options are "mu", or "sigma".
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
fitted.gddlm <- function(object, 
                         type = "mu", 
                         ...) {
    
    if(!inherits(object, "gddlm"))
        stop("model must be of class \"gddlm\".")
    
    type <- match.arg(type, choices = c("mu", "sigma"))
    
    if(type %in% c("mu")) {
        res <- object$mu_fitted
        }
    if(type %in% c("sigma")) {
        res <- object$sigma_fitted
        }
    
    return(res)
    }


#' @rdname fitted.gddlm
#' @export
fitted.values.gddlm <- fitted.gddlm