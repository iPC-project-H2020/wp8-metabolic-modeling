library(caret)
library(tidyverse)
library(MASS)
library(ggpubr)
library(factoextra)

perform_lda_and_silhouette <- function(class.labels, # classification 
                                       data.Matrix, # FBA flux matrix
                                       info.fl, # directory to store silhouette scores
                                       col.sign  #color for the different classes
                                       ){
  
  # Fit the model
  model <- lda(class.labels ~., data = data.Matrix, )
  proportion = model$svd^2/sum(model$svd^2)
  # Make predictions
  predictions <- model %>% predict(data.Matrix)

  dist.lda.space = as.matrix(dist(predictions$x))
  sil = cluster::silhouette(as.numeric(as.factor(data.Matrix$class.labels)), dist.lda.space) 
  sil.width = as.data.frame(sil[1:dim(data.Matrix)[1], ])
  sil.width$category = levels(as.factor(data.Matrix$class.labels))[sil.width$cluster]
  rownames(sil.width) = rownames(data.Matrix)
  write.table(sil.width, 
            file = paste0('./', info_fl, '/silhouette_width_lda.txt'))

  pdf(paste0('./', info_fl, '/lda_silhouette_score.pdf'), width = 10)
  par(mfrow=c(1,2))
  plot(predictions$x, 
     col = col.sign[as.factor(data.Matrix$class.labels)], 
     pch =19)
  #par(xpd=TRUE)
  legend("topright", inset=c(-0.0001,0), 
       col = col.sign[as.factor(data.Matrix$class.labels)],
       legend=unique(data.Matrix$class.labels), 
       pch=19, 
       title="Group", 
       bty = 'o', cex = 0.6)
  plot(sil, main="Silhouette scores")
  dev.off()
  df.description = unlist(strsplit(pr, split = '/'))
  df.description = df.description[-c(1:2)]
  names(df.description) = c('model', 'cl_resolution', 'ic_score_by', 'cl_aggregation_by')
  si.summary = as.numeric(summary(sil)$si.summary)
  names(si.summary) = names(summary(sil)$si.summary)
  silhouette_summary = rbind(silhouette_summary, 
                           c(df.description, 
                             si.summary))
  return(silhouette_summary)
}

spc_cluster_center <- function(class.labels, spca.res){
  mean.cl = lapply(unique(class.labels), function(x){
    sel.lab = which(x==class.labels)
    if (length(sel.lab)==1){
      mean.clust = spca.res$x[names(sel.lab), ]
    } else {
      mean.clust = apply(spca.res$x[names(sel.lab), ], 2, mean)
    }
    return(mean.clust)
  })
  mean.cl = do.call(rbind, mean.cl)
  rownames(mean.cl) <- unique(class.labels)
  return(mean.cl)
}

plot_corr_reactions <- function(spca.res, threshold = NULL, pci, 
                                reaction_subsystem, 
                                react_table = FALSE, 
                                class.labels){
  if (is.null(threshold))   threshold = 1.98*sd(spca.res$vectors[,pci])
  par(mfrow= c(2,2), mar=c(4,6, 2, 2), xpd = NA)
  loadings.pci = spca.res$vectors[order(spca.res$vectors[, pci]), pci]
  loadings.pci.pos = loadings.pci[loadings.pci>threshold]
  loadings.pci.neg = loadings.pci[loadings.pci< -threshold]
  
  col_subSys = colors_signature(reaction_subsystem[names(c(loadings.pci.pos, 
                                                           loadings.pci.neg)),1])
  col.pos = reaction_subsystem[names(loadings.pci.pos), 1]
  col.neg = reaction_subsystem[names(loadings.pci.neg), 1]
  
  barplot(loadings.pci.pos, horiz = TRUE, las=2, 
          main = paste0('sPC ',pci), cex.names = 0.7, cex.axis = 0.7,
          cex.lab = 1, col = col_subSys[col.pos])
  plot(x=1:10, y=1:10, type="n", axes=FALSE, xlab="", ylab="")
  legend("topleft", fill=col_subSys[unique(col.pos)], 
         legend = unique(col.pos), 
         title = "", inset = c(-0.1, 0.1), cex = 1.1, 
         bty = "n")
  
  plot(x=1:10, y=1:10, type="n", axes=FALSE, xlab="", 
       ylab="")
  legend("topleft", fill=col_subSys[unique(col.neg)], 
         legend = unique(col.neg), 
         title = "", inset = c(-0.35, 0.1), cex = 1.1, 
         bty="n")
  barplot(loadings.pci.neg, horiz = TRUE, las=2, 
          main = "", cex.names = 0.7, cex.axis = 0.7,
          cex.lab = 1, col = col_subSys[col.neg])
  if (react_table){
    cl_center = spc_cluster_center(class.labels, spca.res)[, pci]
    gr.pos = cl_center[which.max(cl_center)]
    gr.neg = cl_center[which.min(cl_center)]
    react_table = cbind(c(loadings.pci.pos, loadings.pci.neg), 
                        c(rep(names(gr.pos), length(loadings.pci.pos)),
                          rep(names(gr.neg), length(loadings.pci.neg))))
    return(react_table)
  }else {
    return()
  }
}


## perform supervised PCA with kernel 
perform_spca_and_silhouette <- function(class.labels, # classification 
                                       data.Matrix, # FBA flux matrix
                                       info.fl, # directory to store silhouette scores
                                       col.sign, #color for the different classes
                                       silhouette_summary, 
                                       reaction_subsystem, 
                                       dims = c(1, 2), 
                                       plotSPCA = TRUE 
){
  spca.res = tsiMisc::spca(x = data.Matrix, y = class.labels, center = TRUE, 
                         scale = TRUE, retx = TRUE)
  ncomp = which.min(cumsum(spca.res$values)/sum(spca.res$values)<0.99)
  ## perform kmeans on the reduced space 
  km = kmeans(x = spca.res$x[, 1:ncomp], centers = length(unique(class.labels)), iter.max = 100)

  sign.kmean.cluster = unique(cbind(class.labels[names(km$cluster[order(km$cluster)])], 
                                  km$cluster[order(km$cluster)]))
  xi = dims[1]
  yi = dims[2]
  if (plotSPCA){
    par(mfrow=c(1,1), mar=c(5,4,4,17), xpd=TRUE)
    plot(spca.res$x[, xi], spca.res$x[, yi], pch = 19, col = col.sign[class.labels], 
       xlab = paste0("sPC", xi, " - ",round(spca.res$values[xi]/sum(spca.res$values)*100,2), "%"), 
      ylab = paste0("sPC", yi, " - ",round(spca.res$values[yi]/sum(spca.res$values)*100,2), "%"))
    # points(km$centers[, xi], km$centers[, yi], pch = 15, col = col.sign[sign.kmean.cluster[,1]])
    legend("topright", inset=c(-1.3, 0), legend=names(col.sign), fill = col.sign,
       title="Class", bty = "n", cex = 0.8, xjust = 0, yjust = 0, title.adj = 0.05)
  }
  dist.spca.space = as.matrix(dist(spca.res$x[,1:ncomp]))
  sil = cluster::silhouette(as.numeric(as.factor(class.labels)), dist.spca.space) 
  
  ## store the data 
  data.Matrix$class.labels = class.labels
  data.Matrix$clusters = km$cluster
  sil.width = as.data.frame(sil[1:dim(data.Matrix)[1], ])
  sil.width$category = levels(as.factor(data.Matrix$class.labels))[sil.width$cluster]
  rownames(sil.width) = rownames(data.Matrix)
  write.table(sil.width, 
              file = paste0('./', info_fl, '/silhouette_width_spca.txt'))
  
  pdf(paste0('./', info_fl, '/spca_silhouette_score.pdf'), width = 10)
  par(mfrow=c(1,2))
  plot(spca.res$x[, xi], spca.res$x[, yi], pch = 19, col = col.sign[class.labels], 
       xlab = paste0("sPC", xi, " - ",round(spca.res$values[xi]/sum(spca.res$values)*100,2), "%"), 
       ylab = paste0("sPC", yi, " - ",round(spca.res$values[yi]/sum(spca.res$values)*100,2), "%"))
  #par(xpd=TRUE)
  legend("topright", inset=c(-0.0001,0), 
         col = col.sign[as.factor(data.Matrix$class.labels)],
         legend=unique(data.Matrix$class.labels), 
         pch=19, 
         title="Group", 
         bty = 'o', cex = 0.6)
  plot(sil, main="Silhouette scores")
  dev.off()
  
  pdf(file = paste0('./', info_fl, '/spca_corr_reactions.pdf'))
  react_table = NULL
  for (pci in 1:5){
    react_table_pci = plot_corr_reactions(spca.res, threshold = NULL, pci, 
                        reaction_subsystem, react_table = TRUE, 
                        class.labels = class.labels)
    react_table = rbind(react_table, react_table_pci)
  }
  dev.off()
  
  
  df.description = unlist(strsplit(pr, split = '/'))
  df.description = df.description[-c(1:2)]
  names(df.description) = c('model', 'cl_resolution', 'ic_score_by', 'cl_aggregation_by')
  si.summary = as.numeric(summary(sil)$si.summary)
  names(si.summary) = names(summary(sil)$si.summary)
  silhouette_summary = rbind(silhouette_summary, 
                             c(df.description, 
                               si.summary))
  
  return(list(silhouette_summary=silhouette_summary, react_table = react_table))
}  

myclassification = function(desc.vec){
  attribute = NULL
  if (desc.vec["ic10_plus"]== "high") attribute = paste0(attribute, "_Ews_high")
  if ((desc.vec["ic1_plus"] == "high")|(desc.vec["ic2_plus"] == "high")) attribute = paste0(attribute, "_Proliferating")
  if (desc.vec["ic14_plus"]== "high") attribute = paste0(attribute, "_Glucose/Hypoxia")
  if ((desc.vec["ic6_plus"]== "high")|(desc.vec["ic4_plus"]== "high")) attribute = paste0(attribute, "_Redox/OxPhosPho")
  #if (desc.vec["ic9_plus"]== "high") attribute = paste0(attribute, "_mRNAsplicing")
  if (is.null(attribute)){
    if (desc.vec["ic10_plus"]== "low"){ 
      attribute = paste0(attribute, "_Ews_low")
    } else {
      attribute = "_Average"
    }
  } 
  return(attribute)
}

check.classification <- function(class.labels, pdxScores){
  while(length(which(table(class.labels)==1))>0){
    dist.ic.scores = as.matrix(dist(pdxScores))
    sel.category = names(which(table(class.labels)==1))
    new.cat = lapply(sel.category, function(x){
      sel.pdx.cl = which(class.labels==x)
      nearest.pdx.cl = which(dist.ic.scores[sel.pdx.cl, ] == min(dist.ic.scores[sel.pdx.cl, -sel.pdx.cl ]))
      new.cat <- class.labels[nearest.pdx.cl]
      names(new.cat) = names(sel.pdx.cl)
      return(new.cat)
    })
    new.cat= unlist(new.cat)
    class.labels[names(new.cat)] = new.cat
  }
  return(class.labels)
}

colors_signature = function(class.labels){
  mycols = c(RColorBrewer::brewer.pal(8, "Dark2")[-c(1,5)], 
             RColorBrewer::brewer.pal(10, "Paired"), 
             RColorBrewer::brewer.pal(5, "Set3"))
  mycols = mycols[1:length(unique(class.labels))]
  names(mycols) = unique(class.labels)
  return(mycols)
}

pdx.info.norm.aggregation = function(list.pdx.info){
  grep("PDX", names(list.pdx.info))
  norm.data = lapply(list.pdx.info, function(pdx.element){
    return(pdx.element$pdx.info.norm)
  })
  norm.data = do.call(rbind, norm.data)
  return(norm.data)
}

kmeans.classification = function(info.norm.matrix, ncenters = 5){
  if (length(which(is.na(info.norm.matrix)))!=0) info.norm.matrix[which(is.na(info.norm.matrix))]=0
  pca.info <<- prcomp(info.norm.matrix)
  pdx.cl <<- do.call(rbind, strsplit(rownames(info.norm.matrix), '_'))
  pdx.cl <<- as.data.frame(pdx.cl)
  colnames(pdx.cl) <- c("PDX", 'cluster')
  dist.pca = dist(pca.info$x)
  cl.kmeans <<- kmeans(dist.pca, ncenters)
  names(cl.kmeans$cluster) = rownames(info.norm.matrix)
  p1 <- autoplot(pca.info, loadings = TRUE, loadings.label = TRUE, col=cl.kmeans$cluster)
  p2 <- autoplot(pca.info, data = pdx.cl, col="PDX")
  return(list(cl.kmeans, p1, p2))
}

load_cluster_info <- function(FBA_fn, info_fl){
  cluster_info = lapply(FBA_fn, function(x){
    pdx_cl = unlist(strsplit(x, split= "_"))
    cl.num = unlist(strsplit(pdx_cl[4], split = '.txt'))
    return(c(pdx_cl[2], cl.num))
  })
  cluster_info = do.call(rbind, cluster_info)
  colnames(cluster_info) = c('PDX', 'Cluster')
  cluster_info = as.data.frame(cluster_info)
  
  pdx.info = lapply(unique(cluster_info$PDX), function(x){
    pdx.info = read.csv(paste0('./', info_fl, '/', x, '_cluster_attributes.txt'),
                        sep="")
    mean_val = apply(pdx.info, 2, mean)
    sd_val = apply(pdx.info, 2, sd)
    
    pdx.info.norm = apply(pdx.info, 1, function(y){
      return((y-mean_val)/sd_val)})
    pdx.info.norm = t(pdx.info.norm)[, -dim(pdx.info)[2]]
    rownames(pdx.info.norm) = paste0(x, "_", pdx.info$cluster)
    
    pdx.attribute = apply(pdx.info, 1, function(x){
      sel.high = x > mean_val + 1*sd_val
      sel.low = x < mean_val - 1*sd_val
      sel.avr = (x < mean_val + 1*sd_val)&(x > mean_val - 1*sd_val)
      x[sel.high] = "high"
      x[sel.low] = "low"
      x[sel.avr] ="average"
      return(x)
    })
    return(list(pdx.info = pdx.info,
                pdx.attribute = t(pdx.attribute[-dim(pdx.attribute)[1], ]), 
                pdx.info.norm = pdx.info.norm)
    )
  })
  names(pdx.info) = unique(cluster_info$PDX)
  add.info = apply(cluster_info, 1, function(x){
    pdx = as.character(x[1])
    cl = x[2]
    p.info = pdx.info[[pdx]]$pdx.info
    sel = which(p.info$cluster == cl)
    return(pdx.info[[pdx]]$pdx.attribute[sel,])
  })
  cluster_info = cbind(cluster_info, t(add.info))
  return(list(pdx.info=pdx.info, cluster_info= cluster_info))
}

sd_compare <- function(pdxScores){
  try_class = apply(pdxScores, 2, function(ic_scores){
    aux = ic_scores
    aux[which(ic_scores > sd(ic_scores))] = 'high'
    aux[which(ic_scores < -sd(ic_scores))] = 'low'
    aux[which((ic_scores > -sd(ic_scores)) & 
                (ic_scores < sd(ic_scores)))] = 'average'
    return(aux)
  }) 
  return(try_class)
} 


intersects <- function(rxn_list, out = NULL){
  if (length(rxn_list)<2) return(out)
  len <- seq(length(rxn_list))
  if (missing(out)) out <- list()
  for (idx in split((inds <- combn(length(rxn_list), 2)), col(inds))){
    ii <- len > idx[2] & !(len %in% idx)
    out[[(n<-paste(names(rxn_list[idx]), collapse = "."))]] <- intersect(rxn_list[[idx[1]]], rxn_list[[idx[2]]])
    out <- intersects(append(out[n], rxn_list[ii]), out=out)
  }
  out
}