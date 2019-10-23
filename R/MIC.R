##' @name MIC
##' @rdname MIC
##' @title MIC - Compute a matrix of implied causality
##'
##' @description
##' Generates a Matrix of Implied Causation for a given model, following the
##'   rules in Brick and Bailey (submitted).
##'
##' @param model An SEM model (see details for compatibility)
##' @param latents Compute causal influences of/on latent variables? (default TRUE)
##' @param se Compute SEs? (default TRUE, but only available for some models)
##' @param standardized Compute causal influences using the standardized model parameters? (default FALSE)
##' @param exogenous Compute causal influences of exogenous variables? (default TRUE)
##' @param ... Other parameters passed to mxEval (e.g. defvar.row=)
##'
##' @return The MIC for this model.
##'
##' @details
##' MICs trace all paths from all variables to all others, and present the
##' result as an asymmetric matrix.  The value of a single path is the product
##' of all single-headed arrows taken as part of that path, and the total causal
##' influence is the sum of all paths that connect one variable to another by
##' following only single-headed arrows, and only in the direction of the arrow.
##'
##' Standard errors are computed by the difference method using the mxSE function.
##'
##' Currently accepts MxRAMModel objects, and converts basic lavaan and blavaan objects.
##' Standard Errors are not currently available for converted or standardized models.
##'
##' @seealso MICTable
##'
##' @import OpenMx
##'
##' @export
MIC <- function(model, latents=TRUE, standardized=FALSE, exogenous=TRUE, se=NA, ...) {

    # Convert to MxRAMModel
    model <- as.MxRAMModel(model, standardized=standardized, exogenous=exogenous)

    # Add I and compute:
    model <- mxModel(model, mxMatrix("Iden", nrow(model$A), name="I"))
    solution <- mxEval(solve(I-A), model, ...)
    if(!latents) {
	      Fmat <-mxEval(F, model, ...)
	    solution <- Fmat %&% solution
    }

    outVal <- data.frame(solution)
    ses <- NULL

    # Prep SEs if available:
    if(is.na(se) || se) {

      ses <- tryCatch(expr = suppressWarnings(mxSE(solve(I-A), model, silent=TRUE, ...)),
                      error=function(e) {return(NULL)})
      if(!is.na(se) && is.null(ses)) {
        warning("Could not compute standard errors. If these were not originally OpenMx models, that's why. Look for updates to this in the future, but also see umxLav2RAM() in the meantime.  If they were all MxModels, try running mxCheckIdentification() on your models to make sure they are well-identified.  ")
      }
    }
    if(!is.null(ses)) {
      if(!latents) {
        ses <- Fmat %&% ses
      }

      colnames(ses) <- colnames(solution)
      rownames(ses) <- rownames(solution)
      ses <- data.frame(ses)
   }

   if(!is.null("ses")) {attr(outVal, "SEs") <- ses}
   attr(outVal, "model") <- model$name
   attr(outVal, "MIC") <- TRUE
   return(outVal)
}


##' @name MICor
##' @rdname MICor
##' @title MICor - Matrix of implied correlation for a RAM model in OpenMx.
##'
##' @description
##' Computes a correlation or covariance matrix for a RAM model.  You can do this
##' better by using `mxGetExpected(model, "covariance")`.  This function is here
##' because a reader asked that we make it.
##'
##' @param model An Mx model in the RAM formulation
##' @param latents Include latent variables? (default TRUE)
##' @param stdize Return a correlation matrix instead of a covariance? (default TRUE)
##'
##' @return The covariance or correlation matrix implied by this model.
##'
##' @importFrom stats cov2cor
##'
##' @details
##' See Paper.
computeCorrelation <- function(model, latents=TRUE, stdize=TRUE) {
	  if(!is(model, "MxRAMModel")) {
	      stop("The MICr package currently only analyzes models in the RAM framework.",
	         "If you would like MICs for other frameworks, please provide feedback ",
		          "on the OpenMx forums at https://openmx.ssri.psu.edu/")
	   }
    A <- mxEval(A, model)
    S <- mxEval(S, model)
    I <- diag(1, nrow(A))
    outval <- solve(I-A) %&% S
    if(!latents) {
	      outval <- mxEval(F, model) %&% outval
  }
    if(stdize) { outval <- cov2cor(outval)}
	    return(outval)
}
