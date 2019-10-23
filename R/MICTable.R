##' @name flattenMIC
##' @rdname flattenMIC
##' @title flattenMIC - Create
##'
##' @description
##' Generates a Matrix of Implied Causation for a given model, following the
##'   rules in Brick and Bailey (submitted).
##'
##' @param mic An MIC matrix as returned by \code{\link{MIC}}
##' @param model the (character) label of the model, if different from what's included in the MIC
##' @param cullTiny whether or not to remove very small values from the table (default TRUE)
##' @param tinyCutoff how small is "very small" for culling (default 1e-6)
##' @param includeModel whether or not to include a separate column for model name in the table (default FALSE)
##'
##' @return The MIC table flattened form of this MIC
##'
##' @importFrom methods is
##' @importFrom tidyr pivot_longer
##' @importFrom dplyr select
##' @details
##'  Creates a flat-table version of a single MIC.
flattenMIC <- function(mic, model=NULL, cullTiny=TRUE, tinyCutoff=1e-6,
                       includeModel=FALSE) {

  # Check for pre-flattened MIC
  if(is.data.frame(mic) && "from" %in% names(mic)) {
    return(mic)
  }

  # Get variable names
  mic$to <- setdiff(names(mic), "model")

  # Handle model name if included in MIC
  modelName <- attr(mic, "model")  # NULL if missing
  if(!is.null(model)) {
    modelName <- model
  }
  if(!("model" %in% names(mic))) {
    mic$model <- modelName  # Does nothing if NULL
  }

  # Gather the matrix into a table
  mdf <- mic %>% tidyr::pivot_longer(-one_of("to", "model"), names_to="from",
                                     values_to="value")

  # Remove unwanted entries:
  if(cullTiny) { # If requested, small entries
    mdf <- mdf[abs(mdf$value) > tinyCutoff,]
  }
  mdf <- mdf[mdf$to != mdf$from, ]  # loopback/identical

  # Sort for aesthetic reasons (order keeps ties in order)
  mdf <- mdf[order(mdf$value, decreasing=TRUE),c("model", "from", "to", "value")]
  rownames(mdf) <- NULL

  if(!includeModel) { # If requested, remove the model name
    mdf <- dplyr::select(mdf, -model)
  }

  return(mdf)
}

##' @name MICTable
##' @rdname MICTable
##' @title MICTable
##'
##' @description
##' Takes a set of models and compares them all, following the rules in Brick and Bailey (2019).
##'
##' @param ... MxModels from which to extract MICs
##' @param minAbs the smallest absolute size of causal effect that is worth noticing (default .01)
##' @param minDiff the smallest difference in effect that is worth considering a difference of scale (default minAbs)
##' @param splitByType If TRUE, return a list of tables split by the type of difference (e.g. existence of causal pathway, sign of causal path, scale of path, none of the above).  If FALSE, return a single table. Default TRUE.
##' @param standardize If TRUE, return standardized values
##' @param print A list of the tables to be printed.  Set to FALSE to stop pretty-printing. By default, prints effects that differ by existence, sign, or scale. See details.
##' @param from A list of sources: paths not from a listed source are not printed.  Default NA (print all).
##' @param to A list of the outcomes; paths not ending at a listed outcome are not printed. Default NA (print all).
##' @param se Whether or not to return standard errors (default TRUE).  Model must be a run MxModel, but consult umxLav2RAM to create and run an MxModel using lavaan syntax
##'
##' @return If splitByType is FALSE, a single data frame containing all implied causal effects, as compared across all models.
##'
##' If splitByType is TRUE (the default if more than one model is povided), five data frames:
##'  $exist contains causal effects that exist (|>minAbs|) in at least one model but do not exist in at least one other
##'  $sign contains causal effects that exist in all models, but with different signs in at least two models
##'  $scale contains causal effects with the same sign in all models, but which differ in scales (diff > minDiff) between at least two models
##'  $other causal effects that are identical between models
##'  $all contains the combined table, as would be returned if splitByType were FALSE.
##'
##' @details
##' Returns a comparison table comparing MICs from several models.
##'
##' The print= argument can be used to trigger pretty-printing using kable from the knitr package.
##'   The value of the print= vector contains the list of tables to pretty-print. Each table printed will be reduced to the paths of interest, as defined by the from= and to= arguments.
##'   Possible values are:
##'       existence: Model-implied causal effects with absolute value > minAbs in at least one model that either do not exist or have absolute value < minAbs in at least one other
##'       sign: Model-implied causal effects that exist in all models, but are positive in at least one model and negative in at least one other
##'       scale: Model-implied causal effects that exist with the same sign in all models, but differ in scale by at least minDiff in two models
##'       other: Model-implied causal effects that do not differ by at least minDiff in any model
##'       all: All model-implied causal effects
##'       any: Model implied causal effects that differ that differ by sign, scale, or existence between at least two models
##'    Pretty-printing may be disabled by setting it to FALSE or an empty vector.
##'    If only one model is provided, "all" will print the entire model table; any other legal value will print all effects with an absolute value larger than minAbs.
##' This functionality may be replicated using the kable() function if more precision is desired.
##'
##' @examples
##' XYZmodel <- lavaan("Y ~ .2*X +.6*Z
##'                     Z ~ .8*X")
##' # Print all paths
##' MICTable(XYZmodel, print="all")
##'
##' # Print only sign, scale, and exist as separate tables:
##' MICTable(XYZmodel, print="all")
##'
##' # Print only the nonzero influences of X:
##' MICTable(XYZmodel, from="X")
##'
##' @importFrom rlang list2 enexprs .data
##' @importFrom purrr pmap_lgl
##' @importFrom tidyr pivot_wider
##' @importFrom stats dist
##' @importFrom knitr kable
##' @import dplyr
##' @export
 MICTable <- function(..., minAbs=.01, minDiff=NA, splitByType=TRUE,
                           standardize=FALSE, caption=NULL,
                           print=c("exist", "sign", "scale"),
                           from=NA, to=NA, se=NA) {

  # Handle input
  if(is.na(minDiff)) {minDiff <- minAbs}

  # Processing print argument
  matches <- pmatch(tolower(print), c("existence", "sign", "scale", "all"))

  # Build up model names for table
  modelNames <- as.character(rlang::enexprs(...))
  models <- list2(...)
  for(mNo in seq_along(models)) {
    if(!is.null(attr(models[[mNo]], "model"))) names(models)[mNo] <- attr(models[[mNo]], "model")
    #  Replaces names with MxModel names if they exist.
    else if(is(models[[mNo]], "MxRAMModel")) {
      if(!(models[[mNo]]$name %in% names(models) ||
         startsWith(models[[mNo]]$name, "untitled"))) {
        names(models)[mNo] <- models[[mNo]]$name
      }
    } else { names(models)[mNo] <- as.character(modelNames[[mNo]]) }
  }
  modelNames[!is.na(names(models))] <- names(models)[!is.na(names(models))]

  # Default caption
  if(is.null(caption)) {
    caption <- paste("Implied Causation Table:", paste(modelNames, collapse=", "))
  }

  # Assemble all MICs
  tMICs <- NULL
  tSEs <- NULL
  for(mNo in seq_along(models)) {
    aModel <- models[[mNo]]
    aName <- names(models)[mNo]
    if(is.null(attr(aModel, "MIC"))) {aModel <- MIC(aModel, standardized = standardize, se=se)}
    flatModel <- flattenMIC(aModel, includeModel=TRUE, model=aName, cullTiny=FALSE)
    flatSE <- NULL
    if((is.na(se) || se) && !is.null(attr(aModel, "SE")))
      {flatSE <- flattenMIC(attr(aModel, "SE"), includeModel=TRUE, model=paste0(aName, "_SE"), cullTiny=FALSE)}
    tMICs <- rbind(tMICs, flatModel)
    tSEs <- rbind(tSEs, flatSE)
  }

  se <- (is.na(se) || se) && (!is.null(tSEs) && nrow(tSEs) > 0)  # No SEs if we have no SEs.
  # Widen--one column per model, one row per path
  wideMIC <- pivot_wider(tMICs, names_from="model", values_from="value")
  if(se) {
    wideSE <- pivot_wider(tSEs, names_from="model", values_from="value")
    wideMIC <- left_join(wideMIC, wideSE, by=c("from", "to"))

    # reorder:
    newNameOrder <- intersect(paste(rep(modelNames, each=2), c("", "_SE")), names(wideMIC))
    newNameOrder <- c("from", "to", newNameOrder)
    if(length(newNameOrder) == length(names(wideMIC))) {
      # In case of duplicate naming problems, skip reordering.
      wideMIC <- wideMIC[,newNameOrder]
    }
  }


  # Handle default from and to now that we know what our options are
  if(OpenMx:::single.na(from)) from <- unique(wideMIC$from)
  if(OpenMx:::single.na(to)) to <- unique(wideMIC$to)

  if(length(setdiff(c(from, to), unique(c(wideMIC$from, wideMIC$to)))) > 0) {
    warning(paste("Implied cause or effect of elements",
                  paste(setdiff(c(from, to), unique(c(wideMIC$from, wideMIC$to))), collapse=", "),
                  "was requested, but they do not appear in the MIC.")
            )
  }

  if(length(modelNames) <= 2) {
    # If "all" is requested for pretty-print, print all
    if(5 %in% matches) {
      wideMIC <- wideMIC[wideMIC$from %in% from & wideMIC$to %in% to,]
      print(kable(wideMIC, caption=caption))
      return(invisible(wideMIC))
    } else if(any(!is.na(matches))) {
    # If anything else is requested, pretty print all non-zeros
      wideMIC <- wideMIC[wideMIC$from %in% from & wideMIC$to %in% to & ifelse(is.na(wideMIC[,3]), FALSE, abs(wideMIC[,3] > minAbs)),]
      print(kable(wideMIC, caption=caption))
      return(invisible(wideMIC))
    }
    # If it's set to FALSE, print nothing, and return normally.
    return(wideMIC[wideMIC$from %in% from & wideMIC$to %in% to,])
  }

  # Filter out cases where one variable doesn't exist in a given model
  invalids <- wideMIC %>% dplyr::filter_all(any_vars(is.na(.data)))
  others <- wideMIC %>% dplyr::filter_all(all_vars(!is.na(.data)))

  existenceDifference <- function(...) { length(unique(abs(c(...)) > minAbs)) > 1 }
  signDifference <- function(...) { length(unique(sign(c(...)))) > 1}
  scaleDifference <- function(...) { any(dist(abs(c(...))) > minDiff)}

  # Differ by existence of an effect
  differs <- pmap_lgl(others[,modelNames,drop=FALSE], existenceDifference)
  existence <- others %>% dplyr::filter(differs)
  #%>%                arrange(desc(abs(max(diff(.data[[,modelNames]])))))
  # Do not differ by existence
  others <- others %>% dplyr::filter(!differs)

  # Differ by sign of the effect:
  differs <- pmap_lgl(others[,modelNames,drop=FALSE], signDifference)
  sign <- others %>% dplyr::filter(differs)
  # %>% arrange(desc(max(abs(diff(.data[,modelNames])))))

  # Neither of those:
  others  <- others %>% dplyr::filter(!differs)

  # Differ by scale of the effect:
  differs <- pmap_lgl(others[,modelNames,drop=FALSE], scaleDifference)
  scale <- others %>% dplyr::filter(differs)
  #%>% arrange(desc(max(abs(diff(.data[[,modelNames]])))))


  # Essentially the same:
  others <- others %>% dplyr::filter(!differs) %>% data.frame()
  #%>% arrange(desc(max(abs(diff(.data[[,modelNames]])))))

  # Combinations
  any <- cbind(rbind(existence, sign, scale), type=c(rep("exist", nrow(existence)), rep("sign", nrow(sign)), rep("scale", nrow(scale))))
  all <- rbind(any, cbind(others, type=rep("None", nrow(others))))

  output <- all
  if(splitByType)
    output <- list(exist=existence, sign=sign, scale=scale, other=others, all=all, any=any)
  theTable <- NULL
  if(1 %in% matches) { theTable <- rbind(theTable, existence[existence$from %in% from & existence$to %in% to,]); caption=paste(caption, "Existence &", sep=" ")}
  if(2 %in% matches) { theTable <- rbind(theTable, sign     [sign     $from %in% from & sign     $to %in% to,]); caption=paste(caption, "Sign &", sep=" ")}
  if(3 %in% matches) { theTable <- rbind(theTable, scale    [scale    $from %in% from & scale    $to %in% to,]); caption=paste(caption, "Scale &", sep=" ")}
  if(4 %in% matches) { theTable <- rbind(theTable, other    [other    $from %in% from & other    $to %in% to,]); caption=paste(caption, "Existence &", sep=" ")}
  if(5 %in% matches) { theTable <- rbind(theTable, all      [all      $from %in% from & all      $to %in% to,]); caption=paste(caption, "All &", sep=" ")}
  if(length(matches) > 0) {caption <- substr(caption, 1, nchar(caption) -1)}
  if(nrow(theTable) > 0) {print(kable(theTable, caption=caption))}

  if(any(!is.na(matches))) return(invisible(output))
  else return(output)

}

