% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paretoFront.R
\name{paretoFront}
\alias{paretoFront}
\title{paretoFront - plot the tradeoff of cost and effect}
\usage{
paretoFront(mic, outcome, costs, mapping = NULL, ..., nameCol = "name",
  costCol = "cost", scaleCol = "scale")
}
\arguments{
\item{mic}{An MIC matrix as returned by \code{\link{MIC}}, or a model that can be coerced into one}

\item{outcome}{the (character) label of the outcome}

\item{costs}{A data frame containing columns for the name and cost of each intervention.  Optionally, a column for the scale of the intervention (see details)}

\item{mapping}{additional aesthetics to pass to ggplot (e.g. color)}

\item{...}{additional arguments reserved for later use.}

\item{nameCol}{the column in the data frame containing intervention names (default: "name")}

\item{costCol}{the column in the data frame containing intervention costs (default: "cost")}

\item{scaleCol}{column in the data frame containing intervention scales (default: "scale").  See details.}
}
\value{
A list of Pareto-optimal points, which maximize the cost/effect tradeoff
}
\description{
Plots the pareto front--the best choice at each cost value
}
\details{
Plots the Pareto front, which shows the maximal tradeoff between cost and efficiency of an intervention.

 If the column named in scaleCol (default "scale") is included, this is interpreted to represent the scale of the intervention.
 Simple multiplication will be used to compute the overall effect of this change, assuming linearity in effect.
 This approach makes it possible to examine different elements of intervention (e.g. 3 sessions vs 12 sessions) that may have varying costs.
}
