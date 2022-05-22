bind.cluster.function <- function(enrichment.results, col.names = c("cluster", "function")){
  summary.matrix <- NULL
  for(cl in 1:length(enrichment.results)){
    res.cl = enrichment.results[[cl]]
    if (is.null(res.cl)) res.cl <- matrix(NA, ncol = length(col.names)-1, 
                                          nrow = 10)
    res.cl = cbind(rep(cl-1, 10), res.cl)
    colnames(res.cl) = col.names
    summary.matrix <- rbind(summary.matrix, res.cl)
  }
  colnames(summary.matrix) <- col.names
  return(summary.matrix)
}

pdx_clustering = function(working_dir, input_file,  
                          model_dir, info_dir, 
         clustering_resolution = 0.8, # standard seurat parameter
         determine_data_dim = FALSE, # then 50 dimensions will be taken 
         get.Enrichment = FALSE
         ){
  library(Seurat)
  fn.name = unlist(strsplit(input_file, split = paste0(working_dir, "scRNAseq_data/")))[2]
  fn.name = unlist(strsplit(fn.name, split= ".txt"))
  
  pdx.data <- read.csv(input_file, sep="", row.names = 1)
  
  # Remove rows without hugo identifier
  sel.na <- which(is.na(rownames(pdx.data)))
  if (length(sel.na)!=0) pdx.data <- pdx.data[-sel.na,]
  
  gene.names <- rownames(pdx.data)
  pdx.seurat <- CreateSeuratObject(counts = pdx.data, 
                                   project = fn.name, 
                                   min.cells = 3, 
                                   min.features = 100)
  ## cell filtering
  pdx.seurat <- subset(pdx.seurat, subset = nFeature_RNA > 100)
  ## Find variable genes
  pdx.seurat <- Seurat::FindVariableFeatures(pdx.seurat, selection.method = "vst", nfeatures = 1000)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(pdx.seurat), 10)
  # scale the data before dim reduction 
  all.genes <- rownames(pdx.seurat)
  pdx.seurat <- Seurat::ScaleData(pdx.seurat, features = all.genes)
  # pca
  pdx.seurat <- Seurat::RunPCA(pdx.seurat, features = VariableFeatures(object = pdx.seurat))

  # determine the dimensionality of the data set
  if (determine_data_dim){
    pdx.seurat <- JackStraw(pdx.seurat, num.replicate = 100, dims = 100)
    pdx.seurat <- ScoreJackStraw(pdx.seurat, dims = 1:50)
  }
  # clustering
  pdx.seurat <- FindNeighbors(pdx.seurat, dims = 1:50)
  pdx.seurat <- FindClusters(pdx.seurat, resolution = clustering_resolution)
  
  # UMAP
  pdx.seurat <- RunUMAP(pdx.seurat, dims = 1:50, n.neighbors = 10)
  
  # find markers for every cluster
  pdx.markers <- FindAllMarkers(pdx.seurat, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  tab.markers <- pdx.markers %>%
    group_by(cluster) %>%
    top_n(n = 100, wt = avg_log2FC)
  if (get.Enrichment){ 
    enrichment.results = lapply(unique(tab.markers$cluster), function(cl){
      gostres <- gost(query = tab.markers$gene[tab.markers$cluster==cl], 
                    organism = "hsapiens")
      return(gostres$result$term_name[1:10])
    })
  
    enrichment.results.full = lapply(unique(tab.markers$cluster), function(cl){
    gostres <- gost(query = tab.markers$gene[tab.markers$cluster==cl], 
                    organism = "hsapiens")
    return(gostres$result[1:10, ])
  })
  
    enrichment.results = bind.cluster.function(enrichment.results)
  
    cn = colnames(enrichment.results.full[[1]])
    enrichment.results.full = bind.cluster.function(enrichment.results.full, 
                                                  col.names = c("cluster", cn))
  
    write.table(enrichment.results, file=paste0(info_dir, fn.name, "_cluster_enrichment.txt"), 
              row.names = FALSE)
    write.table(enrichment.results, file=paste0(info_dir, fn.name, "_cluster_enrichment_full.txt"), 
              row.names = FALSE)
  }
  write.table(tab.markers, file=paste0(info_dir, fn.name, "_cluster_markers.txt"), 
              row.names = FALSE)
 
  # save seurat object 
  saveRDS(pdx.seurat, file = paste0(info_dir, fn.name, ".rds"))
  return(pdx.seurat)
}
  
cluster_profile_and_classification = function(working_dir, data.rds, 
                                              gene_number, fl.name, dest_dir, 
                                              score_by, 
                                              aggregation_by) {
  
  ### load signatures ### 
  signature_files_name = list.files(path = paste0(working_dir, "IC_genes/"))
  ic_genes = lapply(signature_files_name, function(fn){
    return(readxl::read_excel(paste0(working_dir, "IC_genes/", fn)))
  })
  ic_component  = unlist(strsplit(signature_files_name, split='.xlsx'))
  names(ic_genes) = ic_component
  
  # compute cluster profiles
  
  score_method = get(score_by)
  aggregation_method = get(aggregation_by)
  
  data_gene_info = remove_genes_not_in_CBM(data.rds, gene_number)
  norm.data = data_gene_info[[1]]
  data_gene_number = data_gene_info[[2]]
  
  cl.profile.per.cluster <- NULL
  cluster_attributes <- NULL
  for (cl in levels(data.rds$seurat_clusters)){
    sel.cl <- which(cl == data.rds$seurat_clusters)
    data.cl <- norm.data[, sel.cl]
    cl.profile <- apply(data.cl, 1, aggregation_method)
    
    idx_component = lapply(ic_genes, function(icg){
      return(data_gene_number[icg$GENE])
    })
    
    ic_scores = lapply(idx_component, function(idx){
      return(score_method(cl.profile[as.character(idx[!is.na(idx)])]))
    })
    
    write.table(cl.profile,
                file = paste0(dest_dir,'/', fl.name, "_cluster_", 
                              cl, ".txt")) 
    cl.profile.per.cluster = cbind(cl.profile.per.cluster, cl.profile)
    cluster_attributes = rbind(cluster_attributes, unlist(ic_scores))
  }
  
  cluster_number = levels(data.rds$seurat_clusters)
  cluster_attributes = cbind(cluster_attributes, cluster = cluster_number)
  
  return(list(pdx_cluster_profile_matrix = cl.profile.per.cluster , 
               cluster_attributes = cluster_attributes))
}

remove_genes_not_in_CBM = function(data.rds, gene_number){
  norm.data <- data.rds@assays$RNA
  data_gene_number <- gene_number[rownames(norm.data)]
  # remove genes not in the CBM and change gene id
  sel.rm <- which(is.na(data_gene_number))
  norm.data <- norm.data[-sel.rm, ]
  rownames(norm.data) <- data_gene_number[-sel.rm]
  data_gene_number <- data_gene_number[-sel.rm]
  return(list(norm.data, data_gene_number))
}


pdx_profiles = function(dest_dir, working_dir, 
                        input_file,  
                        model_dir, 
                        model_name, 
                        clustering_resolution = 0.8, # standard seurat value
                        score_by = "mean",
                        aggregation_by = "median", 
                        determine_data_dim = FALSE, # then 50 dimensions will be taken 
                        getEnrichment = FALSE
){
  
  dest_dir = paste0(dest_dir, "/cl_resolution_", clustering_resolution, 
                    "/ic_score_by_", score_by,
                    "/cl_aggregation_by_", aggregation_by)
  fl.name = unlist(strsplit(input_file, split = paste0(working_dir, "scRNAseq_data/")))[2]
  if (!grepl('.txt', fl.name)) input_file = paste0(input_file, ".txt")
  
  if (length(list.files(path = dest_dir, pattern = fl.name))!=0){
    return(cat('Data already pre-processed for ', fl.name, 
               ' with clustering resolution ', clustering_resolution, 
               ' and model ', model_name,' \n'))
  }
  if (!file.exists(dest_dir)){
    dir.create(dest_dir, recursive= TRUE)
  }
  info_dir = paste0(working_dir, "cluster_info/pdx/", model_name, 
                    "/cl_resolution_", clustering_resolution, 
                    "/ic_score_by_", score_by,
                    "/cl_aggregation_by_", aggregation_by, 
                    "/")
  if (!file.exists(info_dir)){
    dir.create(info_dir, recursive = TRUE)
  }
  
  ## gene symbol to gene GEM identifier
  model.gene.description <- read.table(file = paste0(model_dir, "/gene.tsv"), sep="\t", header = TRUE)
  gene_number <- model.gene.description$gene_number
  
  sel.na = which(model.gene.description$symbol == "")
  if (length(sel.na)!=0){
    gene_number = gene_number[-sel.na]
    gene_symbol = model.gene.description$symbol[-sel.na]
  } else {
    gene_symbol = model.gene.description$symbol
  }
  names(gene_number) <- gene_symbol
  
  ### perform clustering 
  data.rds = pdx_clustering(working_dir, input_file,  
                 model_dir, 
                 info_dir, 
                 clustering_resolution = clustering_resolution, 
                 determine_data_dim = determine_data_dim, 
                 get.Enrichment = getEnrichment)
  
  ######## remove genes not in the CBM #####
  cl.data = remove_genes_not_in_CBM(data.rds, gene_number)
  norm.data = cl.data[[1]]
  data_gene_number = cl.data[[2]]
  
  ################## CLUSTER CLASSIFICATION HERE ##############################
  profile.class = cluster_profile_and_classification(working_dir, data.rds, 
                                     gene_number, fl.name, dest_dir, 
                                     score_by, 
                                     aggregation_by)

  write.table(profile.class$cluster_attributes,
              file = paste0(info_dir, fl.name,
                            "_cluster_attributes.txt")) 
  
  trans_matrix <- profile.class$pdx_cluster_profile_matrix
  colnames(trans_matrix) <- paste0(fl.name, "_", levels(data.rds$seurat_clusters))
  
  trans_matrix <- cbind(data_gene_number,
                        names(data_gene_number),
                        trans_matrix)
  count_matrix <- cbind(data_gene_number,
                        names(data_gene_number),
                        as.matrix(norm.data))
  
  colnames(count_matrix)[1:2] <- colnames(trans_matrix)[1:2] <- c("gene_number", 
                                                                  "symbol")
  write.table(trans_matrix, file = paste0(info_dir, 
                                          fl.name, "_clusters.txt"), row.names = FALSE) 
  write.table(count_matrix, file = paste0(info_dir, 
                                          fl.name, "_all.txt"), row.names = FALSE)
  
}
                        