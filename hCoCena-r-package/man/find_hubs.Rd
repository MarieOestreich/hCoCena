% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_hubs.R
\name{find_hubs}
\alias{find_hubs}
\title{Find Hub Genes}
\usage{
find_hubs(
  clusters = c("all"),
  top = 10,
  tree_layout = F,
  TF_only = F,
  save = F,
  plot = F
)
}
\arguments{
\item{clusters}{Either "all" (default) or a vector of cluster colours for which the hub detection should be performed.}

\item{top}{An integer. All genes are ranked based on their hub potential, this parameter defines the number of the top ranked genes to be considered a hub gene. Default value is 10.}

\item{tree_layout}{A Boolean. Whether or not to depict the network witht tree layout, implying a sort of hierarchical structure to the network.}

\item{TF_only}{Either FALSE (default, all genes in clsuter are considered for hub genes), or "all" (all genes from transcriptionfactor supplementary file are considered for hub genes),
or any gene category listed in the last column of the provided transcriptionfactor supplementary file (only that subgroup condired for hub genes).}

\item{save}{A Boolean. Whether or not a labelled hub-network per cluster and the expression heatmap are to be save to PDF. Default is FALSE.}

\item{plot}{A Boolean. Wheather or not to plot the network (per cluster) with highlighted hub nodes. Default is FALSE.}
}
\description{
Hub genes are determined per cluster using a combined ranking based on weighted degree centrality, weighted closeness centrality and weighted betweenness centrality.
A table of hub genes per cluster is returned as an excel file and a heatmap of the hub genes expression values is plotted per cluster and per dataset.
}
