#' Initialise the Working Directory
#' 
#' init_wd sets up from which directory to read the count data, the annotation data, the reference files and it defines in which directory to create a save folder that all outputs will be written to.
#' @param dir_count_data A string specifying the path of the folder containing the count data files. All count files should be found in the same directory. For further specifications regarding the count files, see the documentation of read_data(). 
#' NOTE: if you are NOT loading your expression data from files, but instead use dataframe already existing in your environment, set this parameter to FALSE.
#' @param dir_annotation A string specifying the path of the folder containing the annotation data files. All annotation files should be found in the same directory. For further specifications regarding the annotation files, see the documentation of read_data(). 
#' NOTE: if you are NOT loading your expression data from files, but instead use dataframe already existing in your environment, set this parameter to FALSE.
#' @param dir_reference_files A string specifying the path of the folder containing the reference files needed for some analysis steps. The files can be found on the GitHub repository under 'reference files'. For details on the files and how to exchange them, refer to the README in the repository.
#' @param dir_output A string specifying the path of the folder in which hCoCena will create a new save folder for the analysis outputs.
#' @seealso [read_data()]
#' @export


init_wd <- function(dir_count_data, dir_annotation, dir_reference_files, dir_output){
	
	hcobject[["working_directory"]][["dir_count_data"]] <<-  dir_count_data
	hcobject[["working_directory"]][["dir_annotation"]] <<-  dir_annotation
	hcobject[["working_directory"]][["dir_reference_files"]] <<-  dir_reference_files
	hcobject[["working_directory"]][["dir_output"]] <<-  dir_output
	
}