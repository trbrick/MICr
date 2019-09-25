##' @name paretoFront
##' @rdname paretoFront
##' @title paretoFront - plot the tradeoff of cost and effect
##'
##' @description
##' Plots the pareto front--the best choice at each cost value
##'
##' @param mic An MIC matrix as returned by \code{\link{MIC}}, or a model that can be coerced into one
##' @param outcome the (character) label of the outcome
##' @param costs A data frame containing columns for the name and cost of each intervention.  Optionally, a column for the scale of the intervention (see details)
##' @param mapping additional aesthetics to pass to ggplot (e.g. color)
##' @param ... additional arguments reserved for later use.
##' @param nameCol the column in the data frame containing intervention names (default: "name")
##' @param costCol the column in the data frame containing intervention costs (default: "cost")
##' @param scaleCol column in the data frame containing intervention scales (default: "scale").  See details.
##'
##' @return A list of Pareto-optimal points, which maximize the cost/effect tradeoff
##'
##' @details
##'  Plots the Pareto front, which shows the maximal tradeoff between cost and efficiency of an intervention.
##'
##'  If the column named in scaleCol (default "scale") is included, this is interpreted to represent the scale of the intervention.
##'  Simple multiplication will be used to compute the overall effect of this change, assuming linearity in effect.
##'  This approach makes it possible to examine different elements of intervention (e.g. 3 sessions vs 12 sessions) that may have varying costs.
##'
##' @importFrom ggplot2 ggplot geom_point geom_text geom_step scale_x_continuous ggtitle aes aes_string expand_scale
##'
##' @export
##'
paretoFront <- function(mic, outcome, costs, mapping=NULL, ...,
                        nameCol="name", costCol="cost", scaleCol="scale") {

  # Handle models and MICs
  if(is.null(attr(mic, "MIC"))) {
    mic <- MIC(mic)
  }
  attr(mic, "model") <- "effect"

  # Make the MICTable with the appropriate elements
  costs[,nameCol] <- as.character(costs[,nameCol])
  ixns <- unique(costs[,nameCol])
  mt <- MICTable(mic, print=FALSE, from=ixns, to=outcome)

  # Total costs for each intervention
  totalCost = merge(costs, mt, by.x=nameCol, by.y="from")

  if(".DisplayName" %in% names(costs)) {
    stop(".DisplayName is reserved by paretoFront. Sorry.")
  }

  totalCost[,".DisplayName"] <- totalCost[,nameCol]

  # Handle optional scaling
  if(scaleCol %in% names(totalCost)) {
    scale <- totalCost[,scaleCol]
    scaled <- scale != 1
    totalCost[scaled, ".DisplayName"] <- paste(totalCost[scaled, nameCol], scale[scaled], sep="@")
    totalCost[scaled, "effect"] <- totalCost[scaled, "effect"] * scale[scaled]
  }

  # Compute Pareto front: sort by cost with max effect first. Keep a cumulative max effect--the first to each max is the front
  # The code to compute the Pareto front is inspired by https://stackoverflow.com/a/21297086.
  # Credit goes to anonymous user6291 and to BrodieG as an indirect effect
  frontCost = totalCost[order(totalCost[,costCol], -totalCost$effect,decreasing=FALSE),]
  frontCost$optimal <- !duplicated(cummax(frontCost$effect))

  frontCost <- frontCost[,c(".DisplayName", nameCol, costCol, "effect", "optimal")]

  # Helper for the final lines:
  front <- frontCost[frontCost$optimal==TRUE,]
  front <- rbind(list("", "", min(front[,costCol]), 0, TRUE),
                 front,
                 list("", "", max(frontCost[,costCol])+.1*diff(range(frontCost[,costCol])), max(front$effect), TRUE))

  # Plot using ggplot.  Inspired by https://stackoverflow.com/q/21296795 and user6291.
  p <- ggplot(frontCost, aes_string(x=costCol, y="effect", color="optimal")) +
    geom_point() +
    geom_text(aes(label=.DisplayName), show.legend=FALSE) +
    geom_step(data=front, direction = "hv") +
    scale_x_continuous(expand=expand_scale(mult=.1)) +
    ggtitle("Pareto Front Showing Optimal Interventions At Each Cost")
  print(p)

  return(frontCost)

}

