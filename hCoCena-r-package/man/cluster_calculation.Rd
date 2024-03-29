% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_calculation.R
\name{cluster_calculation}
\alias{cluster_calculation}
\title{Cluster Calculation}
\usage{
cluster_calculation(
  cluster_algo = "cluster_leiden",
  no_of_iterations = 2,
  resolution = 0.1,
  partition_type = "ModularityVertexPartition",
  max_cluster_count_per_gene = 1,
  return_result = F
)
}
\arguments{
\item{cluster_algo}{The clustering algorithm to be used. The choice is between "cluster_leiden" (default), "cluster_louvain", "cluster_label_prop", "cluster_fast_greedy",
"cluster_infomap", "cluster_walktrap" and "auto" (in which case all are tested and the one with the highest modularity is chosen).}

\item{no_of_iterations}{Some of the algorithms are iterative (e.g. Leiden Algorithm). Set here, how many iterations should be performed.
For information on which other algorithms are iterative, please refer to their documentation in the igraph or leidenbase package. Default is 2.
@param max_cluster_count_per_gene The maximum number of different clusters a gene is allowed to be associated with during the different iterations before it is marked as indecisive and removed.
Default is 1.}

\item{resolution}{The cluster resolution if the cluster algorithm is set to "cluster_leiden". Default is 0.1. Higher values result in more clusters and vice versa.}

\item{partition_type}{Name of the partition type. Select from 'CPMVertexPartition', 'ModularityVertexPartition', 'RBConfigurationVertexPartition' and 'RBERVertexPartition' Default is 'ModularityVertexPartition.}
}
\description{
The function offers several community detection algorithms to identify dense regions in the co-expression network.
These dense regions represent collections of highly co-expressed genes likely to form a functional group.
}
