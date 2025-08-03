#' @title Print a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function prints out some details of the fitted `gddlm` object in a human-readable format. This includes the tuning parameter used, the formula, estimated coefficients, and the gamma divergence.
#'
#' @param x A fitted `gddlm` object.
#' @param digits The number of significant digits to print. Defaults to `max(3, getOption("digits") - 3)`.
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
#' @md
print.gddlm <- function(x,
                        digits = max(3, getOption("digits") - 3),
                        ...) {

    if(!inherits(x, "gddlm"))
        stop("model must be of class \"gddlm\".")

    cat("Double linear model fitted using gamma divergence with gamma =",
        x$gamma_tuning,
        "\n")

    cat("\nCall:",
        paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n",
        sep = "")

    pretty <- function(x) {
        print.default(format(x, digits = digits),
                      print.gap = 2,
                      quote = FALSE)
        }

    if (length(coef.gddlm(x, type = "mu"))) {
        cat("Coefficients for mu (identity link):\n")
        pretty(coef.gddlm(x, type = "mu"))
        }
    else {
        cat("No coefficients for mu.\n")
        }

    cat("\n")

    if(length(coef.gddlm(x, type = "sigma"))) {
        cat("Coefficients for sigma (log link):\n")
        pretty(coef.gddlm(x, type = "sigma"))
        }
    else {
        cat("No coefficients for sigma.\n")
        }

    cat("\n")
    cat("Gamma divergence: ", as.numeric(x$gamma_divergence))
    cat("\n")

    invisible(x)
    }
