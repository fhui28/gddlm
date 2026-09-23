#' @title Extract values from robust model selection for double linear models using gamma divergence.
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Extracts estimates of coefficients characterizing the mean (mu) and the log of the standard deviation (sigma) of the double linear model from a
#' [gddlm_lsa()] fit.
#'
#' @param object An object of class `gddlm_lsa`.
#' @param lambda The tuning parameter at which estimates are estimated. This can either be a string taking values e.g., `AIC`, `BIC`, `EBIC` and `EBIC2`, corresponding to taking the estimates based on minimizing of these information criterion (see [gddlm_lsa()] for more information on these criterion), or a numeric value corresponding to specific tuning parameter value. Similar to [glmnet::coef.glmnet()], in the case of the latter it uses linear interpolation to construct parameter estimates for values of `lambda` that do not coincide with those used in the fitting algorithm.
#' @param ... Not used.
#'
#'
#' @details
#' See associated manuscript (currently in preparation).
#'
#'
#' @return A list with the following elements:
#' \item{mu_coefficients}{A vector of estimated coefficients corresponding to the mean (mu) component of the robust double linear model. Due to the use of the BAR penalty, some of the estimated coefficients will be shrunk to zero.}
#' \item{sigma_coefficients}{A vector of estimated coefficients corresponding to the log standard deviation (sigma) component of the robust double linear model. Due to the use of the BAR penalty, some of the estimated coefficients will be shrunk to zero.}
#' \item{lambda}{The numerical value of the tuning parameter corresponding to the extracted estimates. If `object$control$sequential == TRUE`, then separate `mu_lambda` and `sigma_lambda` are returned.}
#'
#'
#' @author Francis K.C. Hui <francis.hui@anu.edu.au>
#'
#'
#' @examples
#' \donttest{
#' # See help function for the gddlm_lsa() function for an example
#' }
#'
#'
#' @aliases estimates estimate.coque
#' @export
#' @md

estimates <- function(object, ...) {
    UseMethod("estimates")
    }

#' @rdname estimates
#' @export
estimates.gddlm_lsa <- function(object, lambda, ...) {
    if(!inherits(object, "gddlm_lsa"))
        stop("object must be of class 'gddlm_lsa'.")

    ##-------------------
    ## Extract coefficients
    ##-------------------
    if(object$control$sequential == FALSE) {
        if(!is.numeric(lambda)) {
            lambda <- match.arg(lambda, choices = c("AIC", "BIC", "EBIC", "EBIC2"))
            sel_lambda <- which.min(object$information[,lambda])

            sel_mu_coefficients <- object$mu_coefficients_path[,sel_lambda]
            sel_sigma_coefficients <- object$sigma_coefficients_path[,sel_lambda]
            final_lambda <- object$lambda[sel_lambda]
            }

        if(is.numeric(lambda)) {
            lambda_list <- .lambda_interp(lambda_seq = object$lambda, s = lambda)

            sel_mu_coefficients <- object$mu_coefficients_path[, lambda_list$left, drop = FALSE] %*% Diagonal(x = lambda_list$frac) + object$mu_coefficients_path[, lambda_list$right, drop = FALSE] %*% Diagonal(x = 1 - lambda_list$frac)
            sel_mu_coefficients <- as.vector(sel_mu_coefficients)
            names(sel_mu_coefficients) <- rownames(object$mu_coefficients_path)

            sel_sigma_coefficients <- object$sigma_coefficients_path[, lambda_list$left, drop = FALSE] %*% Diagonal(x = lambda_list$frac) + object$sigma_coefficients_path[, lambda_list$right, drop = FALSE] %*% Diagonal(x = 1 - lambda_list$frac)
            sel_sigma_coefficients <- as.vector(sel_sigma_coefficients)
            names(sel_sigma_coefficients) <- rownames(object$sigma_coefficients_path)

            final_lambda <- lambda
            }
        }

    if(object$control$sequential == TRUE) {
        if(!is.numeric(lambda)) {
            lambda <- match.arg(lambda, choices = c("AIC", "BIC", "EBIC", "EBIC2"))
            sel_mu_lambda <- which.min(object$mu_information[,lambda])
            sel_sigma_lambda <- which.min(object$sigma_information[,lambda])

            sel_mu_coefficients <- object$mu_coefficients_path[,sel_mu_lambda]
            sel_sigma_coefficients <- object$sigma_coefficients_path[,sel_sigma_lambda]
            final_mu_lambda <- object$mu_lambda[sel_mu_lambda]
            final_sigma_lambda <- object$sigma_lambda[sel_sigma_lambda]
            }
        if(is.numeric(lambda)) {
            mu_lambda_list <- .lambda_interp(lambda_seq = object$mu_lambda, s = lambda)
            sigma_lambda_list <- .lambda_interp(lambda_seq = object$sigma_lambda, s = lambda)

            sel_mu_coefficients <- object$mu_coefficients_path[, mu_lambda_list$left, drop = FALSE] %*% Diagonal(x = mu_lambda_list$frac) + object$mu_coefficients_path[, mu_lambda_list$right, drop = FALSE] %*% Diagonal(x = 1 - mu_lambda_list$frac)
            sel_mu_coefficients <- as.vector(sel_mu_coefficients)
            names(sel_mu_coefficients) <- rownames(object$mu_coefficients_path)

            sel_sigma_coefficients <- object$sigma_coefficients_path[, sigma_lambda_list$left, drop = FALSE] %*% Diagonal(x = sigma_lambda_list$frac) + object$sigma_coefficients_path[, sigma_lambda_list$right, drop = FALSE] %*% Diagonal(x = 1 - sigma_lambda_list$frac)
            sel_sigma_coefficients <- as.vector(sel_sigma_coefficients)
            names(sel_sigma_coefficients) <- rownames(object$sigma_coefficients_path)

            final_mu_lambda <- lambda
            final_sigma_lambda <- lambda
            }
        }


    ##-------------------
    ## Format
    ##-------------------
    out <- list(mu_coefficients = sel_mu_coefficients,
                sigma_coefficients = sel_sigma_coefficients)
    if(object$control$sequential == FALSE)
        out$lambda <- final_lambda
    if(object$control$sequential == TRUE) {
        out$mu_lambda <- final_mu_lambda
        out$sigma_lambda <- final_sigma_lambda
        }
    return(out)
    }




#' @title A hidden function to use linear interpolation to construct parameter estimates for values of the tuning parameter that do not coincide with those used in the fitting algorithm.
#'
#' @param lambda_seq is the index sequence that is produced by the model. This *must* be a strictly decreasing sequence.
#' @param s is the new vector at which evaluations are required.
#'
#' @return A list containing the left and right indices to take, and the corresponding fraction. The new parameter values are interpolated between the two using the fraction.
#'
#' @note Modified from and acknowledgments goes to the [glmnet::glmnet()] team!
#' @note Note: lambda_seq decreases. you take: \code{sfrac*left+(1-sfrac*right)}
#' @noRd
#' @noMd

.lambda_interp <- function(lambda_seq, s) {

    if(length(lambda_seq) == 1){
        nums <- length(s)
        left <- rep(1,nums)
        right <- left
        sfrac <-rep(1,nums)
    }

    else {
        k <- length(lambda_seq)
        sfrac <- (lambda_seq[1] - s)/(lambda_seq[1] - lambda_seq[k])
        lambda_seq <- (lambda_seq[1] - lambda_seq)/(lambda_seq[1] - lambda_seq[k])

        sfrac[sfrac < min(lambda_seq)] <- min(lambda_seq)
        sfrac[sfrac > max(lambda_seq)] <- max(lambda_seq)

        coord <- stats::approx(lambda_seq, seq(lambda_seq), sfrac)$y
        left <- floor(coord)
        right <- ceiling(coord)

        sfrac <- (sfrac - lambda_seq[right])/(lambda_seq[left] - lambda_seq[right])
        sfrac[left == right] <- 1
        sfrac[abs(lambda_seq[left] - lambda_seq[right]) < .Machine$double.eps] <- 1
        }

    return(list(left = left,
                right = right,frac = sfrac))
    }

