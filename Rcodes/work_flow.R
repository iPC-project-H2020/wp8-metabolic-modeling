library(Seurat)
working_dir = "C:/Users/mpetrizz/OneDrive - INSTITUT CURIE/Metabolism/OrganizedPipeline/github/"
source(paste0(working_dir, "Rcodes/packages.R"))

job_data <- read_excel(paste0(working_dir, "/job_profiles.xlsx"))

for (ii in 1:dim(job_data)[1]){
  x = job_data[ii, ]
  x = as.matrix(x)
  file.name = x[,"file.name"]
  is.cell_line = eval(parse(text=x[,"is.cell_line"]))
  model_name = x[,"model_name"]
  clustering_resolution = as.numeric(x[,"clustering_resolution"])
  input_type_clustering = x[,"input_type_clustering"]
  ic_score_by = x[,"ic_score_by"]
  aggregation_by = x[,"aggregation_by"]

  input_file = paste0(working_dir, "scRNAseq_data/", file.name)
  model_dir = paste0(working_dir, "GEM/", model_name) 

  if(is.cell_line){ 
  
   if (!(input_type_clustering %in% c('ic10', 'temporal'))){
    cat("Error")
    break
   }

    dest_dir = paste0(working_dir, "cluster_profiles/cell_line/", model_name, "/", input_type_clustering)
    source(paste0(working_dir, 'Rcodes/cell_line_profiles_per_', input_type_clustering, '_signature.R'))
    cell_line_profiles(working_dir, input_file, model_dir, dest_dir, model_name,
                                ic_score_by = ic_score_by, 
                                aggregation_by = aggregation_by)

  } else {
   dest_dir = paste0(working_dir, "cluster_profiles/pdx/", model_name)
   source(paste0(working_dir, "Rcodes/pdx_clustering.R"))
   pdx_profiles(dest_dir, working_dir, 
               input_file,  
               model_dir, 
               model_name,
               clustering_resolution = clustering_resolution, # standard seurat value
               score_by = ic_score_by,
               aggregation_by =  aggregation_by, 
               determine_data_dim = FALSE, getEnrichment = FALSE)
  }
}

setwd(working_dir)
system("matlab -nodisplay -r \" run('./MATLABcodes/workflow.m'); exit\"")




