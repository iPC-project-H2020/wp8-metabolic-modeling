### Compute sPCA and silhouette score
rm(list=ls()) # empty the workspace
setwd(' ') ##set the wd to the parent folder
source(file = './Rcodes/functions.R') ## read the R file with functions
model.name = ' ' ## select the initial model

## set data folders
cp_dir = list.dirs("cluster_profiles")
FBA_dir = list.dirs("FBA_solution")
vtr_dir = list.dirs("ValToReact")
info_dir = list.dirs("cluster_info")

cp_dir = cp_dir[grepl('cl_aggregation_by_', cp_dir)]
FBA_dir = FBA_dir[grepl('cl_aggregation_by_', FBA_dir)]
vtr_dir = vtr_dir[grepl('cl_aggregation_by_', vtr_dir)]
info_dir = info_dir[grepl('cl_aggregation_by_', info_dir)]

reaction_subsystem <- read.csv2(paste0("./GEM/", model.name, "/reaction_subsystem.txt"), header=FALSE)
colnames(reaction_subsystem) <- c("Subsystem", "Rxns")
rownames(reaction_subsystem) <- reaction_subsystem$Rxns

cp_dir = cp_dir[grepl(model.name,cp_dir)]

silhouette_summary = NULL
react_table = NULL

for (pr in cp_dir){
  if (!(grepl('pdx', pr))) next
  desc = unlist(strsplit(pr, split='pdx/'))
  FBA_fl = FBA_dir[grep(desc[2], FBA_dir)]
  vtr_fl = vtr_dir[grep(desc[2], vtr_dir)]
  info_fl = info_dir[grep(desc[2], info_dir)]
  
  FBA_fn = list.files( paste0('./', FBA_fl))
  vtr_fn = list.files( paste0('./', vtr_fl))
  
  FBA_profiles = lapply(FBA_fn, function(x){
    FBA_profile = read.csv(paste0('./', FBA_fl, '/', x), 
             row.names=1)
    return(FBA_profile)
  })
  
  pdx.info = load_cluster_info(FBA_fn, info_fl)
  cluster_info = pdx.info$cluster_info
  pdx.info = pdx.info$pdx.info
  pdxScores = lapply(pdx.info, function(x) return(x[[3]]))
  pdxScores = do.call(rbind,pdxScores)
  pdxScores = as.data.frame(pdxScores)
  
  pdxName = unlist(lapply(strsplit(rownames(pdxScores), split="\\."), function(x) x[1]))
  pdxCol = colors_signature(pdxName)
  try_classification = sd_compare(pdxScores)
  class.labels = apply(try_classification, 1, myclassification)
  class.labels = check.classification(class.labels, pdxScores)
  col.sign = colors_signature(class.labels)
  
  reaction.names = lapply(FBA_profiles, row.names)
  reaction.names = unique(do.call(c, reaction.names))
  FBA_profiles_ordered = lapply(FBA_profiles, function(x){
    fl.matrix = matrix(0, nrow = length(reaction.names), ncol=1)
    rownames(fl.matrix) = reaction.names
    fl.matrix[rownames(x), ] <- as.matrix(x)
    return(fl.matrix)
  })
  FBA_profile_matrix = do.call(cbind, FBA_profiles_ordered)
  colnames(FBA_profile_matrix) = paste(cluster_info$PDX, cluster_info$Cluster, 
                                       sep='_')

  data.Matrix = t(FBA_profile_matrix)
  data.Matrix = as.data.frame(data.Matrix)
  
  # remove constant fluxes
  v0 = apply(data.Matrix, 2, var) 
  data.Matrix = data.Matrix[, -which(v0==0)] 
  
  silhouette_summary = perform_spca_and_silhouette(class.labels, # classification 
                                                    data.Matrix, # FBA flux matrix
                                                    info.fl, # directory to store silhouette scores
                                                    col.sign, 
                                                    silhouette_summary, 
                                                    reaction_subsystem, 
                                                    plotSPCA = FALSE)
  
  react_table = rbind(react_table, silhouette_summary$react_table)
  silhouette_summary = silhouette_summary$silhouette_summary
  
}

silhouette_summary = as.data.frame(silhouette_summary)
silhouette_summary$Median = as.numeric(silhouette_summary$Median)
silhouette_summary$Mean = as.numeric(silhouette_summary$Mean)
silhouette_summary$Max. = as.numeric(silhouette_summary$Max.)

colnames(react_table) <- c("loadings", "Annotation")

write.table(x= silhouette_summary, 
            file = paste0("./cluster_info/pdx/", model.name, 
                          "/silhouette_summary.txt"), 
            row.names = FALSE)
write.table(x= react_table, 
            file = paste0("./cluster_info/pdx/", model.name, 
                          "/clusters_and_correlated_reactions.txt"), 
            row.names = TRUE)

