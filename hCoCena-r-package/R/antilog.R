#' Antilog
#' 
#' Reverses a log operation.
#' @param lx The log value.
#' @param base The base to which the log was performed
#' @export

antilog <- function(lx , base) {
  
  lbx <- lx/base::log(base::exp(1) , base = base)
  result <- base::exp(lbx)
  return(result)
}