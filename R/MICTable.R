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
##'
##' @return Three data frames:
##'  $exist contains causal effects that exist (|>minAbs|) in at least one model but do not exist in at least one other
##'  $sign contains causal effects that exist in all models, but with different signs in at least two models
##'  $scale contains causal effects with the same sign in all models, but which differ in scales (diff > minDiff) between at least two models
##'  $other causal effects that are identical between models
##'
##' @details
##' Returns a comparison table comparing MICs from several models.
##'
##' @importFrom rlang list2 enexprs .data
##' @importFrom purrr pmap_lgl
##' @importFrom tidyr pivot_wider
##' @importFrom stats dist
##' @import dplyr
##' @export
 MICTable <- function(..., minAbs=.01, minDiff=NA, splitByType=FALSE) {

  # Handle input
  if(is.na(minDiff)) {minDiff <- minAbs}

  modelNames <- rlang::enexprs(...)
  models <- list2(...)
  for(mNo in seq_along(models)) {
    if(is(models[[mNo]], "MxModel") &&
        !is(models[[mNo]], "MxModel")) {
      stop(paste0("MICr currently only supports RAM models, but model ",
                  models[[mNo]]$name, " is not a RAM model."))
    }

    # TODO: Figure out what this was for.
    if(is(models[[mNo]], "MxRAMModel")) {
      if(models[[mNo]]$name %in% names(models) ||
         startsWith(models[[mNo]]$name, "untitled")) {
          names(models)[mNo] <- as.character(modelNames[[mNo]])
      } else {
        names(models)[mNo] <- models[[mNo]]$name
      }
    }
  }
  modelNames <- names(models)

  # Assemble all MICs
  tMICs <- NULL
  for(mNo in seq_along(models)) {
    aModel <- models[[mNo]]
    aName <- names(models)[mNo]
    if(is(aModel, "MxRAMModel")) {aModel <- MIC(aModel)}
    flatModel <- flattenMIC(aModel, includeModel=TRUE, model=aName, cullTiny=FALSE)
    tMICs <- rbind(tMICs, flatModel)
  }

  # Widen--one column per model, one row per path
  wideMIC <- pivot_wider(tMICs, names_from="model", values_from="value")

  if(length(modelNames) < 2) {
    if(!splitByType) {
      return(wideMIC)
    } else {
      return(list(all=wideMIC))
    }
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
  others <- others %>% dplyr::filter(!differs)
  #%>% arrange(desc(max(abs(diff(.data[[,modelNames]])))))

  # Combinations
  any <- cbind(rbind(existence, sign, scale), type=c(rep("exist", nrow(existence)), rep("sign", nrow(sign)), rep("scale", nrow(scale))))
  all <- rbind(any, cbind(others, type=rep("None", nrow(others))))

  if(!splitByType) return(all)
  return(list(exist=existence, sign=sign, scale=scale, other=others, all=all, any=any))

}

