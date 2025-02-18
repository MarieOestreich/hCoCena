% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highlight_geneset.R
\name{highlight_geneset}
\alias{highlight_geneset}
\title{Highlight Geneset}
\usage{
highlight_geneset(
  gene_set,
  name = NULL,
  col = "black",
  label_offset = 3,
  plot = T,
  save = T
)
}
\arguments{
\item{gene_set}{A vector of gene symbols (as strings) to be highlighted in the network.}

\item{col}{A string specifying the color with which the nodes of the genes are framed to highlight them. Default is "black".}

\item{label_offset}{A numeric specifying the distance between the labeled genes and the network. Default is 3.}

\item{plot}{A Boolean. Wheather or not to plot the network (per cluster) with highlighted hub nodes. Default is FALSE.}

\item{save}{A Boolean. Whether or not a labelled hub-network per cluster and the expression heatmap are to be save to PDF. Default is FALSE.}
}
\description{
Receives and highlights a list of genes in the integrated network.
}
