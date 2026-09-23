#' @title Summary for a fitted gddlm model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Takes a fitted `gddlm` object and produces some useful summaries of the fitted model. Note the summary uses [mgcv::summary.gam()] in the background to produce the summary, so please see the help file for that function for more details.
#'
#' @param object A fitted `gddlm` object.
#' @param digits The number of significant digits to print.
#' @param signif.stars If `TRUE`, `p`-values are printed with stars indicating significance.
#' @param ... Currently not used.
#'
#' @return A summary of the fitted `gddlm` object.
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
#' @importFrom Matrix bdiag
#' @importFrom mgcv gam gaulss summary.gam
#' @importFrom stats printCoefmat
#' @method summary gddlm
#' @md

summary.gddlm <- function(object,
                          digits = max(3, getOption("digits") - 3),
                          signif.stars = getOption("show.signif.stars"),
                          ...) {

    if(!inherits(object, "gddlm"))
        stop("model must be of class \"gddlm\".")

    posthoc_mgcv_fit_fn <- function(object) {
        dofit <- mgcv::gam(list(object$mu_formula,
                                object$sigma_formula),
                           data = object$model,
                           family = gaulss(),
                           method = "ML")
        return(dofit)
        }
    posthoc_mgcv_fit <- posthoc_mgcv_fit_fn(object = object)
    posthoc_mgcv_fit$coefficients <- c(object$mu_coefficients, object$sigma_coefficients)
    names(posthoc_mgcv_fit$coefficients)[1:length(object$mu_coefficients)] <- paste0("mu_", names(object$mu_coefficients))
    names(posthoc_mgcv_fit$coefficients)[-(1:length(object$mu_coefficients))] <- paste0("sigma_", names(object$sigma_coefficients))
    posthoc_mgcv_fit$Vp <- object$V_full #as.matrix(Matrix::bdiag(object$V_mu_coefficients, object$V_sigma_coefficients))
    posthoc_mgcv_fit$Ve <- object$V_full #as.matrix(Matrix::bdiag(object$V_mu_coefficients, object$V_sigma_coefficients))
    posthoc_mgcv_fit$Vc <- object$V_full #as.matrix(Matrix::bdiag(object$V_mu_coefficients, object$V_sigma_coefficients))

    cat("Double linear model fitted using gamma divergence with gamma =",
        object$gamma_tuning,
        "\n")
    .print2(x = summary(posthoc_mgcv_fit),
            gamma_divergence = object$gamma_divergence,
            digits = digits,
            signif.stars = signif.stars)
    }



#' @noRd
#' @noMd
.print2 <- function(x,
                    gamma_divergence,
                    digits,
                    signif.stars,
                    ...) {
    cat("Formula:\n")

    if (is.list(x$formula))
        for (i in 1:length(x$formula))
            print(x$formula[[i]])
    else
        print(x$formula)

    if (length(x$p.coeff) > 0) {
        cat("\nParametric coefficients (mu and sigma):\n")
        stats::printCoefmat(x$p.table,
                            digits = digits,
                            signif.stars = signif.stars,
                            na.print = "NA",
                            ...)
    }

    cat("\n")
    if(x$m > 0) {
        cat("Significance of smooth terms:\n")
        stats::printCoefmat(x$s.table,
                            digits = digits,
                            signif.stars = signif.stars,
                            has.Pvalue = TRUE,
                            na.print = "NA",
                            cs.ind = 1,
                            ...)
    }

    cat("\n")
    cat("Residual degrees of freedom:", x$residual.df, "\n")
    cat("Gamma divergence:", round(as.numeric(gamma_divergence), digits))
    cat("\n")

    #if(!is.null(x$rank) && x$rank< x$np) cat("Rank: ",x$rank,"/",x$np,"\n",sep="")
    #if(!is.null(x$r.sq)) cat("R-sq.(adj) = ",formatC(x$r.sq,digits=3,width=5),"  ")
    #if (length(x$dev.expl)>0) cat("Deviance explained = ",formatC(x$dev.expl*100,digits=3,width=4),"%",sep="")

    #cat("Scale est. = ",formatC(x$scale,digits=5,width=8,flag="-"),"  n = ",x$n,"\n",sep="")
    invisible(x)
}

