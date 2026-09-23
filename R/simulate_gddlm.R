#' @title Simulate responses from a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Simulate one or more responses from a fitted `gddlm` object.
#'
#' @param object A fitted `gddlm` object.
#' @param nsim Number of responses to simulate. Defaults to 1.
#' @param seed An optional seed for the random number generator. Please see the corresponding argument in [stats::simulate()] for more details.
#' @param ... Currently not used.
#'
#' @return A data frame of simulated responses
#'
#' @note
#' Acknowledgments to Hannes Riebl and the `lmls` package on whose code this function is adapted off.
#'
#' @author
#' Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#' @examples
#' \donttest{
#' # See the main `gddlm` help file for examples.
#' }
#'
#' @importFrom stats rnorm runif simulate
#' @export
#' @md
simulate.gddlm <- function(object,
                           nsim = 1,
                           seed = NULL,
                           ...) {

    if(!inherits(object, "gddlm"))
        stop("model must be of class \"gddlm\".")

    # RNG handling taken from simulate.lm():
    # https://github.com/wch/r-source/blob/master/src/library/stats/R/lm.R

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        runif(1)  # initialize the RNG if necessary
        }

    if (is.null(seed)) {
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        }
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }

    n <- nrow(object$model) * nsim
    out <- stats::rnorm(n, mean = object$mu_fitted, sd = object$sigma_fitted)
    out <- matrix(out, nrow = nrow(object$model), ncol = nsim)

    out <- as.data.frame(out)
    names(out) <- paste("sim", seq_len(nsim), sep = "_")
    attr(out, "seed") <- RNGstate

    return(out)
    }
