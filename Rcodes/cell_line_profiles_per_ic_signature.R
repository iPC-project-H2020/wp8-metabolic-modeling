gene.expression.profile.per.category <- function(sel.cat, 
                                                 aggregation_meth = aggregation_meth, 
                                                 cell.info,
                                                 dat){
  cel.info.per.category <- cell.info[sel.cat, ]
  xx <- unique(cel.info.per.category[, 'type'])
  
  sel.cells <- cel.info.per.category[which(cel.info.per.category[,'type']==xx), 'name']
  df.per.cat <- dat[,sel.cells]
  gene.names <-  rownames(df.per.cat)
  
  final.ds <- matrix(0, nrow=length(unique(gene.names)), ncol=dim(cel.info.per.category)[1])
  rownames(final.ds) <- unique(gene.names)
  colnames(final.ds) <- cel.info.per.category[, 'name']
  final.ds[rownames(df.per.cat), colnames(df.per.cat)] <- as.matrix(df.per.cat)
  cat.profile <- apply(final.ds, 1, aggregation_meth)
  return(cat.profile)
}

geneId_to_GemId <- function(profile, file.name, dest_dir, gene_number){
  ## remove genes that do not have a gene id mapping in the GEM
  data_gene_number <- gene_number[names(profile)]
  sel.rm <- which(is.na(data_gene_number))
  profile <- profile[-sel.rm]
  # change gene id coding
  names(profile) <- data_gene_number[-sel.rm]
  # save ! 
  write.table(profile,
              file = paste0(dest_dir,
                            file.name, ".txt"))
  return(profile)
}

cell_line_profiles = function(working_dir, input_file, model_dir, dest_dir, model_name, 
                                       ic_score_by = "mean", 
                                       aggregation_by = "median", 
                                       plot_density_ic10_score = FALSE){
  dest_dir = paste0(dest_dir, "/ic_score_by_", ic_score_by, "/cl_aggregation_by_", aggregation_by, '/')
  info_dir = paste0(working_dir, "cluster_info/cell_line/", model_name, 
                    "/ic/ic_score_by_", ic_score_by, 
                    "/cl_aggregation_by_", aggregation_by, '/')
  
  if (!file.exists(dest_dir)){
    dir.create(dest_dir, recursive= TRUE)
  }
  if (!file.exists(info_dir)){
    dir.create(info_dir, recursive= TRUE)
  }
  ## sample annotation
  sample_annotation <- read.delim(paste0(working_dir,"/scRNAseq_data/sample_annotation.txt"))
  cell.names = sample_annotation$CELL[sample_annotation$GROUP == "ASP14_TS"]
  
  ## gene symbol to gene GEM identifier
  model.gene.description <- read.table(file = paste0(model_dir, "/gene.tsv"), sep="\t", header = TRUE)
  gene_number <- model.gene.description$gene_number
  names(gene_number) <- model.gene.description$symbol
  ## I10 associated gene set 
  I10_plus_genes <- readxl::read_excel(paste0(working_dir, "IC_genes/ic.xlsx"))
  
  ic_score_action = get(ic_score_by)
  aggregation_meth = get(aggregation_by)
  
  # read the input data
  dat <- read.delim(input_file)
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  dat <- dat[, cell.names]
  
  ## Consider the genes that correlate with I10+
  var_i10_genes <- apply(dat[I10_plus_genes$GENE,], 1, var)
  not.var = which(var_i10_genes==0)
  var.i10.genes = I10_plus_genes$GENE
  if (length(not.var)!=0) var.i10.genes = I10_plus_genes$GENE[-not.var] 
  
  ic10_score = apply(dat[I10_plus_genes$GENE,], 2, ic_score_action, na.rm=TRUE)
  cell.info <- cbind(rep(NA,length(names(ic10_score ))), names(ic10_score ))
  colnames(cell.info) <- c('type', 'name')
  cell.info <- as.data.frame(cell.info)
  
  if (plot_density_ic10_score){
    X11()
    plot(density(ic10_score))}
  
  low.thr <- quantile(ic10_score)[2]
  mid.low.thr <- quantile(ic10_score)[3]
  mid.high.thr  <- quantile(ic10_score)[4]
  
  sel.low <- which(ic10_score < low.thr)
  sel.mid.low <- which((ic10_score >= low.thr)&(ic10_score<mid.low.thr))
  sel.mid.high <-  which((ic10_score >= mid.low.thr)&(ic10_score<mid.high.thr))
  sel.high <-  which(ic10_score >= mid.high.thr)
  
  cell.info$type[sel.low] <- 'low'
  cell.info$type[sel.mid.low] <- 'midlow'
  cell.info$type[sel.mid.high] <- 'midhigh'
  cell.info$type[sel.high] <- 'high'
  
  profile.low <- gene.expression.profile.per.category(sel.low, aggregation_meth, cell.info, dat )
  profile.mid.low <- gene.expression.profile.per.category(sel.mid.low, aggregation_meth, cell.info, dat)
  profile.mid.high <- gene.expression.profile.per.category(sel.mid.high, aggregation_meth, cell.info, dat)
  profile.high <- gene.expression.profile.per.category(sel.high, aggregation_meth, cell.info, dat)
  
  plow = geneId_to_GemId(profile.low, file.name = 'ASP14_i10_plus_low', 
                         dest_dir, gene_number )
  pmlow = geneId_to_GemId(profile.mid.low, file.name = 'ASP14_i10_plus_midlow',
                         dest_dir, gene_number )
  pmhigh= geneId_to_GemId(profile.mid.high, file.name = 'ASP14_i10_plus_midhigh', 
                         dest_dir, gene_number )
  phigh = geneId_to_GemId(profile.high, file.name = 'ASP14_i10_plus_high', 
                          dest_dir, gene_number )
  
  p.matrix = cbind(plow, pmlow, pmhigh, phigh)
  write.table(p.matrix, file =paste0(info_dir, "cell_line_i10.txt"))
}






