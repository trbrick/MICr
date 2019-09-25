##' @name standardizeMxRAMModel
##' @rdname standardizeMxRAMModel
##' @title standardizeMxRAMModel: Create a new, standardized MxRAMModel from an old one.
##'
##' @description
##' Transforms an MxRAMModel into a Standardized MxRAMModel
##'
##' @param model an MxRAMModel
##'
##' @return An MxRAMModel containing the same path structure as the original model, but standardized scores
##'
##' @details
##' This function is experimental, and may have bugs.
##' At present, it does not handle constraints, groups, or pretty much anything else
##'   that's at all complicated.
##'
##' @import OpenMx
##'
standardizeMxRAMModel <- function(model) {

  stdVals <- mxStandardizeRAMPaths(model, SE=FALSE)
  for(idx in 1:nrow(stdVals)) {
    aRow <- stdVals[idx,]
    if(aRow$matrix == "A") {
      model$A$values[aRow$row, aRow$col] <- aRow$Std.Val
    } else if(aRow$matrix == "S") {
      model$S$values[aRow$row, aRow$col] <- aRow$Std.Val
    }
  }

  return(model)

}

##' @name extractModelStructure
##' @rdname extractModelStructure
##' @title Extract Model Structure and create an MxRAMModel
##'
##' @description
##' Extract paths from a model and convert to MxRAMModel
##'
##' @param model a path modeling object (see details for supported types)
##' @param standardized Transform all variables into standardized forms? (default FALSE)
##' @param exogenous Include exogenous variables? (default TRUE)
##'
##' @return An MxRAMModel containing the same path structure as the original model
##'
##' @details
##' This function is experimental, and may have bugs.
##'
##' At present, it does not handle constraints, groups, or pretty much anything else
##'   that's at all complicated.
##'
##' Really, this is just a helper to handle the fact that blavaan doesn't export its blavaan class.
##'
##' Currently supported: OpenMx RAM models (easy!), lavaan and blavaan models
##'
##'
extractModelStructure <- function(model, exogenous=TRUE, standardized=FALSE, ...) {

    if(class(model) == "blavaan") {class(model) <- "lavaan"}

  return(as.MxRAMModel(model, exogenous=exogenous, standardized=standardized, ...))
}

##' @name as.MxRAMModel
##' @rdname as.MxRAMModel
##' @title as.MxRAMModel: Create an MxRAMModel from a lavaan or OpenMx model object
##'
##' @description
##' Transforms a model into an MxRAMModel
##'
##' @param model a path modeling object (see details for supported types)
##' @param standardized Transform all variables into standardized forms? (default FALSE)
##' @param exogenous Include exogenous variables? (default TRUE)
##'
##' @return An MxRAMModel containing the same path structure as the original model
##'
##' @details
##' This function is experimental, and may have bugs.
##' At present, it does not handle constraints, groups, or pretty much anything else
##'   that's at all complicated.
##'
##' Currently supported: OpenMx RAM models (easy!), lavaan and blavaan models
##'
as.MxRAMModel <- function(model, exogenous=TRUE, standardized=FALSE, ...) {

  stop(paste0("Not sure how to extract paths from a model of class ", class(model), ".\n",
              "Try building this as an MxModel or a lavaan model."))
}

##' @name as.MxRAMModel.lavaan
##' @rdname as.MxRAMModel.lavaan
##' @title Extract path structure from a lavaan object
##'
##' @description
##' Extracts the paths from a lavaan object and returns them as an MxRAMModel
##'
##' @param model a lavaan or blavaan object
##' @param standardized Compute causal influences for the standardized model? (default FALSE)
##' @param exogenous Compute causal influences of exogenous variables? (default TRUE)
##' @param ... Not currently used.
##'
##' @return An MxRAMModel containing the same path structure as the provided lavaan model
##'
##' @details
##' This function is experimental, and may have bugs.
##' At present, it does not allow for constraints,
##'
##' @importFrom lavaan parTable lavaanNames lavaan
##' @importClassesFrom lavaan lavaan
##' @import OpenMx
##'
as.MxRAMModel.lavaan <- function(model, exogenous=TRUE, standardized=FALSE, ...) {

  params <- lavaan::parTable(model)

  if(length(unique(params$group)) > 1) {
    warning("The MICr package does not currently support multiple group or multilevel models.",
            "Be cautious in interpreting the results.")
  }

  manifests <- lavaan::lavaanNames(model, type="ov")
  latents   <- lavaan::lavaanNames(model, type="lv")

  # Extract paths from a lavaan model.
  # Note: Index here is a hack to get around apply() semantics
  getPaths <- function(index, paramList, exo) {
    parm <- paramList[index,]
    if(!exo && parm$exo == TRUE) {
      return(NULL)
    }

    if(parm$op == "~1") {
      return(mxPath(from="one", to=parm$lhs,
                    values=parm$est, free=(parm$free!=0), arrows=2,
                    labels=ifelse(nchar(parm$label)==0, NA, parm$label)))
    }

    if(parm$op == "=~") {
      return(mxPath(from=parm$lhs, to=parm$rhs,
                    values=parm$est, free=(parm$free!=0), arrows=1,
                    labels=ifelse(nchar(parm$label)==0, NA, parm$label)))
    }

    if(parm$op == "~~") {
      return(mxPath(from=parm$lhs, to=parm$rhs,
                    values=parm$est, free=(parm$free!=0), arrows=2,
                    labels=ifelse(nchar(parm$label)==0, NA, parm$label)))
    }

    if(parm$op == "~") {
      return(mxPath(from=parm$rhs, to=parm$lhs,
                    values=parm$est, free=(parm$free!=0), arrows=1,
                    labels=ifelse(nchar(parm$label)==0, NA, parm$label)))
    }

    # TODO: Constraints are op %in% c("==", "<", ">")

    return(NULL)

  }

  pathModel <- mxModel(manifestVars=manifests, latentVars=latents, type="RAM",
                       sapply(1:nrow(params), getPaths, paramList=params, exo=exogenous))

  if(standardized) pathModel <- standardizeMxRAMModel(pathModel)

  return(pathModel)

}


setMethod("as.MxRAMModel", signature("lavaan" ), as.MxRAMModel.lavaan)
# setMethod("as.MxRAMModel", signature("blavaan"), as.MxRAMModel.lavaan)  # Causes issues.

##' @name as.MxRAMModel.MxModel
##' @rdname as.MxRAMModel.MxModel
##' @title Extract path structure from an OpenMx model
##'
##' @description
##' Extracts the paths from an OpenMx object and returns them as an MxRAMModel
##'   There's really no reason to use this.  Like, why are you even bothering?
##'
##' @param model an OpenMx model
##' @param standardized Compute causal influences for the standardized model? (default FALSE)
##' @param exogenous Compute causal influences of exogenous variables? (default TRUE)
##'
##' @return An MxRAMModel containing the same path structure as the provided lavaan model
##'
##' @details
##' This function is experimental, and may have bugs.
##' At present, it does not allow for constraints,
##'

as.MxRAMModel.MxModel <- function(model, exogenous=TRUE, standardized=FALSE, ...) {

  if(!is(model, "MxRAMModel")) {
    stop("MICr only currently supports MxRAMModels.  Please create your model with type=\"RAM\".")
  }

  if(standardized) model <- standardizeMxRAMModel(model)

  return(model)

}

setMethod("as.MxRAMModel", signature("MxRAMModel"), as.MxRAMModel.MxModel)
setMethod("as.MxRAMModel", signature("MxModel"   ), as.MxRAMModel.MxModel)


##' @name as.MxRAMModel.character
##' @rdname as.MxRAMModel.character
##' @title Extract path structure from a lavaan model string
##'
##' @description
##' Extracts the paths from a lavaan syntax string and returns them as an MxRAMModel
##'
##' @param model a lavaan model string
##' @param standardized Compute causal influences for the standardized model? (default FALSE)
##' @param exogenous Compute causal influences of exogenous variables? (default TRUE)
##' @param ... other arguments to umxLav2RAM()
##'
##' @return An MxRAMModel containing the same path structure as the provided lavaan model
##'
##' @details
##' This function is experimental, and may have bugs.
##' Under the hood, this function uses the umx function umxLav2RAM().  See the help file for that function for additional help.
##'
##' @importFrom umx umxLav2RAM
##'
as.MxRAMModel.character <- function(model, exogenous=TRUE, standardized=FALSE, ...) {

  tryCatch(umxLav2RAM(model, ...), error=function(x) {stop("The only character strings that can be interpreted as models",
                                                        "by this function are model specs using lavaan syntax.")})

  if(standardized) model <- standardizeMxRAMModel(model)

  return(model)

}
setMethod("as.MxRAMModel", "character", as.MxRAMModel.character)

