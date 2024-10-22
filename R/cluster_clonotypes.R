
#' Cluster Clonotypes By Distribution
#'
#' This function plots a UMAP where cells can be colored according to a variable.
#' Multiple clonotypes will also be highlighted.
#'
#' @param seurat_object Seurat object to operate on
#' @param minimum_count Integer specifying minimum number fo clones for a clonotype to be included in analysis
#' @param avg_samples Boolean that determines whether or not frequencies will be determined on a per-sample basis and then averaged across samples of the same subject. Default is 'False', in which case frequencies will be determined by aggregating samples from the same subject.
#' @param exclude_cell_types list of cell_types to exclude from clustering
#' @param exclude_samples list of cell_types to exclude from clustering
#' @param n_families number of clone families to define via hierarchical clustering
#' @param family_colors colors to label clone families. By default ("NULL"), colors provided to the 'plotting_aesthetics' function will be used
#' @param plot Boolean specifying whether or not to plot a clustered heatmap
#' @return an object with two attributes: 'clone_families' (assignment of each clone to a numerical family/cluster) and 'matrix' (dataframe with clonotypes as rows and their propotional distribution across cell_types as columns)
#' @export
cluster_clonotypes <- function(seurat_object,
                               minimum_count=5,
                               avg_samples = FALSE,
                               exclude_cell_types=c(),
                               exclude_samples=c(),
                               n_families=7,
                               family_colors=NULL,
                               plot=TRUE) {

  output <- c()
  
  clonotype_column <- seurat_object@misc$clonotype_column
  cell_type_column <- seurat_object@misc$cell_type_column
  sample_column <- seurat_object@misc$sample_column
  subject_column <- seurat_object@misc$subject_column
  
  meta <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,]
  
  # Add column to describe subject clonotype combination
  meta$subject_clonotype <- paste(meta[[subject_column]], meta[[clonotype_column]], sep="_")
  
  # Exclude some cell types
  meta <- meta[(meta[[cell_type_column]] %in% exclude_cell_types)==FALSE,]
  meta <- meta[(meta[[sample_column]] %in% exclude_samples)==FALSE,]
  
  if (avg_samples == FALSE) {

    
    ## Pheno Count
    pheno_clone_count <- meta %>% 
      group_by(meta["subject_clonotype"], meta[cell_type_column]) %>%
      summarise(pheno_clone_freq=n()) %>%
      group_by(subject_clonotype) %>%
      mutate(total_clone_count=sum(pheno_clone_freq))
    
    # Filter by clone size
    pheno_clone_count <- pheno_clone_count[pheno_clone_count$total_clone_count >= minimum_count,]
    
    # Get frequency of pheno clone
    pheno_clone_count$pheno_clone_freq <- pheno_clone_count$pheno_clone_freq / pheno_clone_count$total_clone_count
    
    # Reshape to wide
    pheno_clone_count <- pivot_wider(pheno_clone_count, id_cols = c("subject_clonotype","total_clone_count"), names_from = cell_type_column, values_from = "pheno_clone_freq")
    
    # Fill na
    pheno_clone_count[is.na(pheno_clone_count)] <- 0
    
    # Drop total_clone_count column
    pheno_clone_count <- pheno_clone_count[ , !(names(pheno_clone_count) == "total_clone_count")]
    
    
  } else {
    
    clone_count <- meta %>% group_by(subject_clonotype) %>%
      summarise(clone_count=n())
    
    # Filter by clone size
    clones_to_include <- clone_count[clone_count$clone_count >= minimum_count, "subject_clonotype"]
    meta <- meta[meta[,"subject_clonotype"] %in% clones_to_include$subject_clonotype,]
    
    ## Add Pheno Count
    pheno_clone_count <- meta %>% group_by(meta["subject_clonotype"], meta[sample_column], meta[cell_type_column]) %>%
      summarise(pheno_clone_count=n())
    
    # Make Template
    subject_clonotypes <- unique(pheno_clone_count$subject_clonotype)
    cell_types <- unique(pheno_clone_count[[cell_type_column]])
    samples <- unique(pheno_clone_count[[sample_column]])
    
    template <- expand.grid(cell_types, subject_clonotypes, samples)
    colnames(template) <- c(cell_type_column, "subject_clonotype", sample_column)
    
    # Combine template and count data
    pheno_clone_count <- left_join(x = template, y = pheno_clone_count, by = c(cell_type_column, sample_column, "subject_clonotype"))
    pheno_clone_count[is.na(pheno_clone_count)] <- 0
    
    # Reshape to wide
    pheno_clone_count <- reshape(pheno_clone_count, idvar = c("subject_clonotype",sample_column), timevar = cell_type_column, direction = "wide")
    
    # Divide by total number of clones
    pheno_clone_count[,-c(1,2)] <- pheno_clone_count[,-c(1,2)] / rowSums(pheno_clone_count[,-c(1,2)])
    
    # Remove rows with NA. These are samples where a clone was not detected at all.
    ## Phenotype frequencies are only averaged across samples in which the clone is observed
    pheno_clone_count <- pheno_clone_count[complete.cases(pheno_clone_count), ]
    
    # Drop sample column 
    pheno_clone_count <- pheno_clone_count[ , !(names(pheno_clone_count) == sample_column)]
    
    # Average frequencies across samples from the same individual
    pheno_clone_count <- pheno_clone_count %>%
      group_by(subject_clonotype) %>%
      summarise(across(everything(), mean))
    
  }
  
  # Format for output
  pheno_clone_count <- as.data.frame(pheno_clone_count)
  row.names(pheno_clone_count) <- pheno_clone_count$subject_clonotype
  pheno_clone_count <- subset(pheno_clone_count, select=-c(subject_clonotype))
  
  colnames(pheno_clone_count) <- gsub("pheno_clone_count.", "", colnames(pheno_clone_count))
  
  ### FOR RANK BASED DISTANCE
  rank_matrix <- pheno_clone_count
  for (rank in 1:ncol(rank_matrix)) {
    
    # Get values to search for on a per clonotype basis
    value_of_next_in_rank <- apply(pheno_clone_count, 1, function(x) sort(x, decreasing = TRUE)[rank])
    
    for ( cell_type in colnames(rank_matrix) ) {
      
      rank_matrix[pheno_clone_count[cell_type]==value_of_next_in_rank,cell_type] <- rank
      
    }
    
  }
  
  # Deal with ties in first rank (random since they are rare)
  no_preference <- apply(rank_matrix, 1, function(x) length(which(x==1)))
  no_preference <- no_preference[no_preference==0]
  
  no_preference <- rank_matrix[names(no_preference), ]
  
  for ( n in 1:nrow(no_preference) ) {
    
    row <- no_preference[n,]
    
    tie <- colnames(row[,row==min(row)])
    
    rank_matrix[row.names(row), sample(tie, 1)] <- 1
  }
  
  # dist
  for (rank_num in 1:ncol(rank_matrix)) {
    
    rank_matrix_zeroed <- rank_matrix
    
    rank_matrix_zeroed[rank_matrix_zeroed!=rank_num] <- 0
    
    dist_to_add <- rank_matrix_zeroed %>%
      dist() %>%
      as.matrix()
    
    dist_to_add[dist_to_add!=0] <- 1/(3^rank_num)
    
    if (rank_num == 1) {
      distance_matrix <- dist_to_add
    } else {
      distance_matrix <- distance_matrix + dist_to_add
    }
    
    
  }
  
  d <- as.dist(distance_matrix)
  
  ##########################
  
  # Hierarchical clustering
  #d <- dist(as.matrix(rank_matrix), method = "manhattan")
  hc <- hclust(d, method="ward.D2")
  
  # Group families
  clone_families <- cutree(hc, k = n_families)
  
  # Get cell_type order
  ordered_clones <- row.names(pheno_clone_count[hc$order,])
  ordered_fams <- unique(clone_families[ordered_clones])
  
  # Renumber clone_families
  clone_families_renumbered <- c()
  for (i in 1:n_families) {
    
    to_add <- clone_families[clone_families==ordered_fams[[i]]]
    to_add[to_add==ordered_fams[[i]]] <- i
    
    clone_families_renumbered <- c(clone_families_renumbered, to_add)
  }
  
  clone_families <- clone_families_renumbered
  
  cell_type_order <- c()
  for (fam in 1:n_families) {
    clones <- names(clone_families[clone_families==fam])
    pheno_means <- colMeans(pheno_clone_count[clones,])
    fam_pheno <- names(pheno_means[pheno_means==max(pheno_means)])
    
    cell_type_order <- c(cell_type_order, fam_pheno)
  }
  
  # need to add back cell types that dont have their own fam
  add_back <- setdiff(colnames(pheno_clone_count), cell_type_order)
  cell_type_order <- c(cell_type_order, add_back) %>% unique()
  
  # Reorder matrix
  pheno_clone_count <- pheno_clone_count[, cell_type_order]
  
  if (plot==TRUE) {
    
    # Colors for plotting aesthetic
    if (is.null(family_colors)) {

      mean_df <- merge(pheno_clone_count, data.frame(clone_families), by.x = 0, by.y = 0)
      mean_df$Row.names <- NULL
      
      mean_df <<- mean_df %>%
        group_by(clone_families) %>%
        summarise(across(everything(), mean))

      color_keys <<- mean_df$clone_families
      mean_df$clone_families <- NULL
      
      max_cell_types <<- colnames(mean_df)[apply(mean_df,1,which.max)]
      
      color_order <<- match(max_cell_types, seurat_object@misc$cell_type_order)
      
      colColors <<- seurat_object@misc$colors[color_order]
      
    # User provided colors
    } else {
      color_dict <- family_colors
      names(color_dict) <- 1:n_families
      
      color_keys <- clone_families[row.names(pheno_clone_count)]
    
      colColors <- color_dict[color_keys]
    }
    
    # Heatmap
    heatmap(t(as.matrix(pheno_clone_count)), scale='none', Rowv=NA, 
              ColSideColors=colColors, cexRow = 2, Colv=as.dendrogram(hc),
              margins = c(15,15)
            )
    
  }
  
  output$clone_families <- clone_families
  output$matrix <- pheno_clone_count
  
  return(output)

}
