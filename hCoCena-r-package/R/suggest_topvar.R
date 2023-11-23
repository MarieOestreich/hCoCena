#' Calculate Expression Variance And Rank
#' 
#' To determine the most variable genes used as input for network calculations the function identifies the inflection point in a curve of the logged variance of ranked genes. 
#' The calculation of the inflection point is a fast approach to identify a first threshold for potentially interesting genes that describe differences in the data while excluding non-variable genes and thereby significantly reducing the calculation time during further analyses.
#' 
#' @noRd

rank_variance <- function(expr){
  # calculate variance:
  var.df <- base::data.frame(variance = base::apply(expr, 1, stats::var))
  var.df[["gene"]] <- base::rownames(expr)
  # determine ranks base on decreasing variance:
  var.df <- var.df[base::order(var.df[["variance"]], decreasing = T),]
  var.df[["rank"]] <- 1:base::nrow(var.df)
  return(var.df)
}

#' Find Loess Span
#' 
#' @noRd

find_span <- function(var.df){

  for(i in base::rev(base::seq(0.5, 1, 0.01))){
    lo <- stats::loess(base::log(var.df$variance)~var.df$rank, span = i)
    xl <- base::seq(base::min(var.df$rank),base::max(var.df$rank), (base::max(var.df$rank) - base::min(var.df$rank))/1000)
    out = stats::predict(lo,xl)
    infl <- base::as.logical(base::diff(base::sign(base::diff(base::diff(out, differences = 1)))))
    if(base::length(infl[infl == T]) > 0){
      break
    }
  }
  return(i)
}

#' Plot Inflection Points:
#' 
#' @noRd

plot_inflections <- function(var.df, i, setname = "my data"){
  lo <- stats::loess(base::log(var.df$variance)~var.df$rank, span = i)
  xl <- base::seq(base::min(var.df$rank),base::max(var.df$rank), (base::max(var.df$rank) - base::min(var.df$rank))/1000)
  out = stats::predict(lo,xl)
  infl <- base::as.logical(base::diff(base::sign(base::diff(base::diff(out, differences = 1)))))
  plot(var.df$rank, log(var.df$variance), type="l", xlab = "Rank", ylab = "log(Variance)", main = setname)
  points(xl[infl ], out[infl ], col="red")
  print(base::paste0(setname, ": Inflection points at the following #genes: "))
  print(base::ceiling(xl[infl ]))
}


#' Suggest Top Most Variant Genes
#' 
#' For each dataset, calculate the inflection points of the logged variance of gene expression values.
#'  Returned values may be used for "topvar" in set_layer_settings().
#' @export


suggest_topvar <- function(){
  for(l in 1:base::length(hcobject[["layers"]])){
    var.df <- rank_variance(hcobject[["data"]][[base::paste0("set", l, "_counts")]])
    i <- find_span(var.df = var.df)
    plot_inflections(var.df = var.df, i = i, setname = hcobject[["layers_names"]][[l]])
  }
}