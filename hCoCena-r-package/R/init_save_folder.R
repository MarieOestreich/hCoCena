#' Creates a Save Folder
#' 
#' A folder with the given name is created in the output directory. All analysis outputs will be saved to this folder.
#' @param name The of the folder to be created.
#' @export

init_save_folder <- function(name){
  
  save_folder <- name
  
  if(!save_folder %in% base::list.dirs(hcobject[["working_directory"]][["dir_output"]])) {
    
    base::dir.create(paste0(hcobject[["working_directory"]][["dir_output"]], save_folder))
    
  }
  
  hcobject[["global_settings"]][["save_folder"]] <<- name
}