# Helper functions for clustering of data in Seurat objects. 

# Dependencies:
#   - Packages (don't load here, but require them to be loaded into the global env 
#               wherever this script is being sourced):
#     - Seurat
#     - tidyverse
#     - glue

###############################################################################
# Rename clusters by size ####
###############################################################################

#' rename_clusters_by_size
#'
#' @description
#' Renames clusters such that the largest cluster is cluster 1 and the smallest cluster is cluster n. 
#'
#' @param object Seurat object
#' @param idents.to.rename Name of the metadata column which contains the idents to be renamed (character string)
#' @param idents.renamed Name for the new metadata column to be created (character string)
#' @param remove.old.names 
#' @param verbose Whether to delete the metadata column with old names (logical)
#'
#' @return Seurat object

rename_clusters_by_size <- function(object,
                                    idents.to.rename,
                                    idents.renamed,
                                    remove.old.names = FALSE,
                                    verbose = TRUE) {
  idents.to.rename.sym <- sym(idents.to.rename)
  idents.renamed.sym <- sym(idents.renamed)
  
  # order clusters by size ====================================================
  clust_order <- object@meta.data %>% 
    group_by(!!idents.to.rename.sym) %>% 
    tally() %>% 
    ungroup() %>% 
    arrange(-n) %>%
    mutate(!!idents.renamed.sym := 1:n()) %>%
    select((!!idents.to.rename.sym), (!!idents.renamed.sym))
  
  # rename clusters by size ===================================================
  object@meta.data[[idents.renamed]] <- clust_order[[idents.renamed]][match(
    x = as.character(object@meta.data[[idents.to.rename]]),
    table = clust_order[[idents.to.rename]]
  )]
  
  object@meta.data[[idents.renamed]] <- factor(
    object@meta.data[[idents.renamed]],
    levels = unique(sort(object@meta.data[[idents.renamed]]))
  )
  
  # set new idents as active idents ===========================================
  object <- SetIdent(object, value = idents.renamed)
  if(remove.old.names) {
    object@meta.data[[idents.to.rename]] <- NULL
  }
  
  # wrap up ===================================================================
  if(verbose) {
    print(glue("Renamed {nrow(clust_order)} clusters by size."))
  }
  
  return(object)
}

###############################################################################
# Merge small clusters ####
###############################################################################

#' merge_small_clusters
#' @description
#' Merges small clusters into their neighbours
#'
#' @param object Seurat object
#' @param assay Which assay to use (e.g. "RNA", "SCT"). Default is DefaultAssay(object).
#' @param min.cells Merge clusters with fewer than min.cells cells (numeric)
#' @param genes Genes to use for calculating correlation between clusters. Default is VariableFeatures(object)
#' @param idents.to.merge Name of metadata column that contains cluster identities to be merged. Character string.
#' @param rename.clusters.by.size Rename clusters by size (1 = largest, n = smallest). Default is TRUE.
#' @param verbose default is TRUE
#'
#' @return Seurat object with merged clusters in metadata and set to active.ident.

merge_small_clusters <- function(object, 
                                 assay = DefaultAssay(object), 
                                 min.cells = 10, 
                                 genes = VariableFeatures(object), 
                                 idents.to.merge, 
                                 rename.clusters.by.size = TRUE,
                                 verbose = TRUE) {
  time_start <- Sys.time()
  
  # merge clusters ============================================================
  object <- SetIdent(object, value = idents.to.merge)
  
  counter <- 0
  
  ## loop to merge smallest cluster into nearest cluster until all clusters > min.cell
  while (min(table(object@active.ident)) <= min.cells) {
    
    smallest_cluster <- attr(which.min(table(object@active.ident)), "names") # which.min returns the first index if more than one
    
    avg.exp <- AverageExpression(object, 
                                 assays = assay, 
                                 features = genes, 
                                 return.seurat = FALSE, 
                                 verbose = FALSE) %>% 
      .[[assay]] %>% 
      log1p() %>% 
      as.data.frame()
    
    correlations <- cor(x = avg.exp[, names(avg.exp) == smallest_cluster], # between smallest cluster and all other clusters
                        y = avg.exp[, names(avg.exp) != smallest_cluster],
                        method = "pearson") %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column(var = "cluster")
    
    nearest_cluster <- correlations[which.max(correlations$V1), "cluster"]
    new_cluster_name <- paste0(smallest_cluster, ".", nearest_cluster)
    
    # create a new vector for merged clusters and fill in new cluster name
    merged_clusters <- object@active.ident %>% as.character()
    merged_clusters[merged_clusters %in% c(smallest_cluster, nearest_cluster)] <- new_cluster_name
    
    # set new cluster names as active.ident
    merged_clusters <- as.factor(merged_clusters)
    object$cluster_merged_min.cell_temp <- merged_clusters
    object <- SetIdent(object, value = "cluster_merged_min.cell_temp")
    
    if (verbose) {
      print(glue("Merged cluster {smallest_cluster} into cluster {nearest_cluster} (correlation: {corr.print}).", 
                 corr.print = round(correlations[correlations$cluster == nearest_cluster, "V1"], 2)))
    }
    
    counter <- counter + 1
  }
  
  # if any clusters were merged: rename clusters ==============================
  if (counter > 0) {
    
    ## message
    if (verbose) {
      print(glue("Merged {counter} cluster(s) with {min.cells} or fewer cells \\
               into their respective nearest neighbours. Now there are {num.clusters} clusters.",
               num.clusters = object$cluster_merged_min.cell_temp %>% unique() %>% length()))
    }
    
    ## rename clusters
    if (rename.clusters.by.size) {
      object <- rename_clusters_by_size(object = object, 
                                        idents.to.rename = "cluster_merged_min.cell_temp", 
                                        idents.renamed   = "cluster_merged_min.cell", 
                                        remove.old.names = TRUE, verbose = TRUE)
    }
    
  }
  
  
  # if no clusters were merged: ===============================================
  if (counter == 0) {
    if (verbose) {
      print(glue("All clusters already have more than {min.cells} cells. \\
                 No clusters were merged."))
    }
  }
  
  # wrap up ===================================================================
  time_end <- Sys.time()
  time_elapsed <- time_end - time_start
  
  if (verbose) {
    print(glue("Finished. \\
               Elapsed time: {round(time_elapsed[[1]], 2)} {attr(time_elapsed, \"units\")}."))
  }
  
  return(object)
}

###############################################################################
# Merge similar clusters ####
###############################################################################

#' merge_similar_clusters
#' @description
#' Merge nearest neighbour clusters that meet the
#'  similarity requirement based on number of DE genes between them.
#'
#' @param object Seurat object
#' @param assay Assay. Default is DefaultAssay(). (character string)
#' @param genes Vector of genes that will be used for calculating average expression, correlation, and DE analysis. (vector)
#' @param idents.to.merge Metadata column with clusters to merge. (character string)
#' @param idents.merged.name Name for the new clusters (character string)
#' @param corr.quantile.cutoff Quantile threshold for correlation values between pairs
#'   of clusters below which clusters should NOT be merged. Default is 0.50. (double, between 0 and 1)
#' @param num.de.genes Number of DE genes between two clusters, above which clusters 
#'   will not be merged. Default is 20. (integer)
#' @param fold.diff Fold difference for a gene to be considered DE. Default is 2. (double)
#' @param p.adj Adjusted p value below which genes will be considered DE. Default is 0.05. (double between 0 and 1)
#' @param rename.clusters.by.size Whether to rename merged clusters by size. Default is TRUE. (logical)
#' @param verbose 
#'
#' @return Seurat object

merge_similar_clusters <- function(object,
                                   assay = DefaultAssay(object),
                                   genes = VariableFeatures(object),
                                   idents.to.merge,
                                   idents.merged.name,
                                   corr.quantile.cutoff = 0.5,
                                   num.de.genes = 20,
                                   fold.diff = 2,
                                   p.adj = 0.05,
                                   rename.clusters.by.size = TRUE,
                                   verbose = TRUE) {
  time_start <- Sys.time()
  
  # merge similar clusters ====================================================
  object <- SetIdent(object, value = idents.to.merge) # set idents
  
  counter <- 1 # to force the first iteration of the while loop
  while (counter > 0) { # as long as some clusters are being merged
    
    counter <- 0
    
    ## calculate averages
    avg.exp <- AverageExpression(object, 
                                 assays = assay, 
                                 features = genes,
                                 return.seurat = FALSE, 
                                 verbose = FALSE) %>% 
      .[[assay]] %>% 
      log1p()
    
    ## measure correlations and sort in decreasing order so most similar pairs are tested first.
    correlations <- cor(avg.exp, avg.exp)
    correlations[lower.tri(correlations, diag = TRUE)] <- NA
    correlations <- as.data.frame(correlations) %>% 
      rownames_to_column(var = "cluster_1")
    correlations <- pivot_longer(correlations, 
                                 cols = 2:ncol(correlations), 
                                 names_to = "cluster_2", 
                                 values_to = "corr", 
                                 values_drop_na = TRUE)
    correlations <- correlations %>% 
      filter(corr > quantile(corr, {{ corr.quantile.cutoff }})) %>% 
      arrange(-corr)
    
    ## merge clusters
    tested <- c() # Write names of merged clusters to this vector. If these clusters appears in a pair again in the loop, skip the pair. 
    clusters.merged <- as.character(object@active.ident)
    
    for (i in 1:nrow(correlations)) {
      
      if ((!correlations$cluster_1[i] %in% tested) & (!correlations$cluster_2[i] %in% tested)) {
        
        de.genes <- FindMarkers(object,
                                ident.1 = correlations$cluster_1[i],
                                ident.2 = correlations$cluster_2[i],
                                assay = "RNA", 
                                min.pct = 0.1, 
                                features = genes, 
                                base = 2,
                                verbose = FALSE) %>% 
          filter(p_val_adj < {{ p.adj }}, 
                 abs(avg_log2FC) > log({{ fold.diff }}, base = 2))
        
        if (nrow(de.genes) < num.de.genes) {
          
          new_cluster <- paste0(correlations$cluster_1[i], ".", correlations$cluster_2[i])
          clusters.merged[clusters.merged == correlations$cluster_1[i]] <- new_cluster
          clusters.merged[clusters.merged == correlations$cluster_2[i]] <- new_cluster
          
          tested <- append(tested, c(correlations$cluster_1[i], correlations$cluster_2[i]))
          
          if (verbose) {
            print(glue("Merged clusters {correlations$cluster_1[i]} and {correlations$cluster_2[i]} \\
                   (DE genes: {nrow(de.genes)}; correlation: {round(correlations$corr[i], 2)})."))
          }
          
          counter <- counter + 1
          
        }
      }
    }
    
    # set merged clusters as idents
    object@meta.data[["cluster_merged_temp"]] <- clusters.merged
    
    object@meta.data[["cluster_merged_temp"]] <- factor(
      object@meta.data[["cluster_merged_temp"]],
      levels = unique(sort(object@meta.data[["cluster_merged_temp"]]))
    )
    
    object <- SetIdent(object, value = "cluster_merged_temp")
    
  }
  
  # if any clusters were merged: ==============================================
  if (!identical(as.character(object@meta.data[["cluster_merged_temp"]]),
                 as.character(object@meta.data[[idents.to.merge]]))) {
    
    ## message
    if (verbose) {
      print(glue("Finished merging nearest-neighbour clusters with fewer than \\
               {num.de.genes} DE genes with a fold difference of {fold.diff}."))
    }
  }
  
  # if no clusters were merged: ===============================================
  if (identical(as.character(object@meta.data[["cluster_merged_temp"]]),
                as.character(object@meta.data[[idents.to.merge]]))) {
    
    ## message
    if (verbose) {
      print(glue("No cluster pairs were found with fewer than \\
               {num.de.genes} DE genes with a fold difference of {fold.diff}. \\
               No clusters were merged."))
    }
  }
  
  # rename clusters ===========================================================
  # if no clusters were merged, the new clusters will identical to input clusters.
  if (rename.clusters.by.size) {
    object <- rename_clusters_by_size(object = object, 
                                      idents.to.rename = "cluster_merged_temp", 
                                      idents.renamed   = idents.merged.name, 
                                      remove.old.names = TRUE, verbose = TRUE)
  }
  
  # wrap up ===================================================================
  time_end <- Sys.time()
  time_elapsed <- time_end - time_start
  
  if (verbose) {
    print(glue("Finished. \\
               Elapsed time: {round(time_elapsed[[1]], 2)} {attr(time_elapsed, \"units\")}."))
  }
  
  return(object)
}

###############################################################################
# Merge clusters (wrapper)
###############################################################################

#' merge_clusters
#' @description
#' Wrapper that runs `merge_small_clusters()` followed by `merge_similar_clusters()`
#'
#' @param object Seurat object
#' @param idents.to.merge meta.data column containing idents to merge. Passed to `merge_small_clusters()`. (character string)
#' @param min.cells Passed to `merge_small_clusters()`. (integer)
#' @param idents.merged.name Name for merged clusters column. (character string)
#' @param ... Passed to `merge_similar_clusters()`
#'
#' @return Seurat object

merge_clusters <- function(object,
                           idents.to.merge,
                           min.cells = 10,
                           idents.merged.name,
                           ...) {
  
  # merge small clusters ======================================================
  print(glue("Merging small clusters ..."))
  object <- merge_small_clusters(object = object, 
                                 min.cells = min.cells, 
                                 idents.to.merge = idents.to.merge)
  
  # merge similar clusters ====================================================
  print(glue("Merging similar clusters ..."))
  # if any clusters were merged in the last step, merge them further by
  # similarity; else merge original clusters by similarity
  if ("cluster_merged_min.cell" %in% names(object@meta.data)) {
    object <- merge_similar_clusters(object = object, 
                                     idents.to.merge = "cluster_merged_min.cell", 
                                     idents.merged.name = idents.merged.name,
                                     ...)
    
    # remove min.cell merged clusters
    object@meta.data[["cluster_merged_min.cell"]] <- NULL
    
  } else {
    object <- merge_similar_clusters(object = object, 
                                     idents.to.merge = idents.to.merge, 
                                     idents.merged.name = idents.merged.name,
                                     ...)
  }
  
  return(object)
}

# Message =====================================================================
message("Loaded helper function `merge_clusters()`, which internally runs:\n",
        "\t1. merge_small_clusters()\n",
        "\t2. merge_similar_clusters()\n",
        "\t3. rename_clusters_by_size()")

# end =========================================================================