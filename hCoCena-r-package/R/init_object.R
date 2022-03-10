#' Initialise hCoCena Object
#' 
#' Creates an object called 'hcobject' in the global envirentment that will be used by hCoCena throughout the analysis.
#' @export

init_object <- function(envo = .GlobalEnv){
	hcobject <- list("working_directory" = list(),
						"data" = list(),
						"supplementary_data" = list(),
						"global_settings" = list(),
						"layer_settings" = list(),
						"layers" = list(),
						"supplement" = list(),
						"layers_names" = list(),
						"layer_specific_outputs" = list(),
						"integrated_output" = list(),
						"cutoff_vec" = NULL,
						"satellite_outputs" = list())

	base::assign("hcobject", hcobject, envir = envo)
}