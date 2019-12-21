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
##' @param se Compute SEs? (default NA; see details)
##' @param standardized Compute causal influences using the standardized model parameters? (default FALSE)
##' @param exogenous Compute causal influences of exogenous variables? (default TRUE)
##' @param N Number of observations for hypothetical standard errors
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
##' Standard errors are computed by the delta method using the mxSE function. Options
##' for SE computation are FALSE (do not compute), "observed" (compute using existing data),
##' "hypothetical" (compute assuming perfect model fit). Any other value (default) selects observed if
##' the model has been run and has not been modified since running, and hypothetical otherwise.
##' Observed computes the SE as a traditional standard error on the fitted model; this is best if
##' the goal is to examine an existing model fit.
##' Hypothetical SEs assume the model is correct. The function generates a covariance and means
##' matrix from the model expectation and fits a new model using those values as data.  It then
##' computes SEs using the delta method from these observed data.  If data exists in the model,
##' the number of data rows specified there is used.  If not, N must be specified.
##'
##' Currently accepts MxRAMModel objects, and converts basic lavaan and blavaan objects.
##' Standard Errors are not currently available for converted or standardized models.
##'
##' @seealso MICTable
##'
##' @import OpenMx
##'
##' @export
MIC <- function(model, latents=TRUE, standardized=FALSE, exogenous=TRUE, se=NA, N=NA, ...) {

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
    seMatch = pmatch(se, c("observed", "hypothetical", "FALSE"), nomatch=-1)
    if(seMatch !=3) {  # If computing SEs at all
      if(seMatch==1) {
        ses <- computeObservedSEs(model=model, quiet=is.na(se), ...)
      } else if(seMatch==2) { # Hypothetical
        ses <- computeHypotheticalSEs(model, quiet=is.na(se), N = N, ...)
      } else if(seMatch<0)  {
        ses <- computeObservedSEs(model=model, quiet=TRUE, ...)
        if(is.null(ses)) {
          ses <- computeHypotheticalSEs(model, quiet=is.na(se), N = N, ...)
        }
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

##' @name single.na
##' @title single.na
##'
##' @description is.na that appropriately handles things that are not single NAs
##'
##' @param a a thing
##'
##' @return Whether or not this is a single NA
single.na <- function (a)
{
  return((length(a) == 1) && (is.list(a) || is.vector(a) ||
                                is.matrix(a)) && (is.na(a) == TRUE))
}


##' @name computeObservedSEs
##' @title Compute observed standard errors from an mxModel
##'
##' @description Simply computes SEs for all total effects using mxSE, throwing appropriate errors for MICr.
##'
##' @param model The model
##' @param quiet If TRUE, do not throw warnings. (default: FALSE)
##' @param ... Additional parameters for mxSE
##'
##' @return computed SEs, or NULL if computation was not possible.
##'
##' @importFrom OpenMx mxSE
##'
computeObservedSEs <- function(model, quiet=FALSE, ...) {
  # Helper to save typing
  warn <- function(warning) {
    if(!quiet) warning(warning, call.=FALSE)
  }

  # If not an MxModel
  if(!is(model, "MxModel")) {
    warn("SEs are currently only available for mxModel objects.  If these were not originally OpenMx models, that may be why. Look for updates to this in the future, but also see umxLav2RAM() in the meantime.")
    return(NULL)
  }

  # If not fitted with these data.
  if(!model@.wasRun || model@.modifiedSinceRun) {
    warn("MIC observed SEs are only available if the model has been successfully fitted and not modified since fitting; use se=\"hypothetical\" to estimate hypothetical SEs.  If you already were trying to use a hypothetical SE, there may have been an error running your model.")
    return(NULL)
  }

  ses <- tryCatch(expr = suppressWarnings(mxSE(solve(I-A), model, silent=TRUE, ...)),
                  error=function(e) {return(NULL)})
  if(is.null(ses)) {
    warn("Could not compute observed standard errors. If they were all MxModels, try running mxCheckIdentification() on your models to make sure they are well-identified.  ")
  }
  return(ses)
}

##' @name computeHypotheticalSEs
##' @title Compute hypothetical standard errors from an mxModel
##'
##' @description Creates standard errors for the model assuming that the
##' model is exactly correct and complete data is available.
##'
##' @param model The model
##' @param N The number of participants modeled.  If NA (the default) and the model has data, the size of the data is used.
##' @param quiet If TRUE, do not throw warnings (default: FALSE)
##' @param ... Additional parameters for mxSE
##'
##' @return computed SEs, or NULL if computation was not possible.
##'
##' @details Computes the model-implied covariance, mean (if available)
##' and thresholds (if available) for the model, fits the model using
##' those as the data set, and computes standard errors for the MIC accordingly.
##'
##' @importFrom OpenMx mxSE
##'
computeHypotheticalSEs <- function(model, N=NA, quiet=TRUE, ...) {
  # Helper to save typing
  warn <- function(warning) {
    if(!quiet) warning(warning, call.=FALSE)
  }

  # If not an MxModel
  if(!is(model, "MxModel")) {
    warn("MIC SEs are currently only available for mxModel objects.  If these were not originally OpenMx models, that may be why. Look for updates to this in the future, but also see umxLav2RAM() in the meantime.")
    return(NULL)
  }

  nobs <- N   # Get specified N
  if(is.na(nobs) && !is.null(model$data)) {
      nobs <- model$data$numObs  # or, failing that, num from data
  }

  # Or warn.
  if(is.null(model$data) && is.na(nobs)) {
    warn("To compute hypothetical SEs without raw or covariance data already in the model, specify a hypothetical population size with N=.")
    return(NULL)
  }

  # Compute expectations
  expectations <- mxGetExpected(model, c("covariance", "means", "thresholds"))
  expCov <- expectations$covariance
  expMeans <- NA
  if(nrow(expectations$means) > 0) {
    expMeans <- expectations$means
  }
  expThresh <- NA
  if(nrow(expectations$thresholds) > 0) {
    expThresh <- expectations$thresholds
  }

  # Re-fit model
  newMod <- NA
  if(single.na(expMeans) && single.na(expThresh)) {
    newMod <- mxModel(model, mxData(expCov, type="cov", numObs = nobs))
  } else if(single.na(expThresh)){
    newMod <- mxModel(model, mxData(expCov, type="cov", means=expMeans, numObs = nobs))
  } else {
    newMod <- mxModel(model, mxData(expCov, type="cov", means=expMeans, thresholds = expThresh, numObs = nobs))
  }

  newMod <- mxRun(newMod, silent = TRUE, suppressWarnings = TRUE)

  # Compute observed SEs
  return(computeObservedSEs(newMod, quiet, ...))
}

##' @name computeCorrelation
##' @rdname computeCorrelation
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
