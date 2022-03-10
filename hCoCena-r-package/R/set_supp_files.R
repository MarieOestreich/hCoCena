#' Set Supplementary Files
#' 
#' Receives the file names of the supplementary files
#' @param Tf A file name. The file must contain a list of at least two columns. One column should be titled with the organism from which the datasets originate and contain gene names, the other column must always be the last column and contain the information if the gene is a transcription factor ("TF"),
#' 	a co-factor ("Co_factor"), a chromatin remodelling protein ("Chromatin_remodeller") or a ribonucleic acid binding protein ("RNBP"). The name of the last column may vary. An exmemplary file for mouse and human is found in the Reference Files Folder in the repository.
#' @param Hall A file name. DESCRIBE ONLY IF NECESSARY.
#' @param Go A file name for a .gmt gene ontology file. You can find one in the reference file folder in the repository.
#' @export

set_supp_files <- function(Tf, Hall, Go){

	hcobject[["supplement"]] <<- base::c(Tf, Hall, Go)

}