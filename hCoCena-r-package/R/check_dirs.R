#' Fixes Directories 
#' 
#' Iteratively calls fix_dir() on the provided working directory paths to fix them if necessary or throw an error if they are invalid.
#' @export

check_dirs <- function(){
  for(x in base::names(hcobject[["working_directory"]])){
    hcobject[["working_directory"]][[x]] <<- fix_dir(hcobject[["working_directory"]][[x]])
  }
}

