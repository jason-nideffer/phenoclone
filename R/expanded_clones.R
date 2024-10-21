
#' Identifying Expanded Clones
#'
#' This function takes a seurat object with meta data that already includes
#' TCR annotations and calculates the expansion of clones in sample.1 compared to sample.2
#'
#' @param seurat_object seurat object to operate on
#' @param sample.1 Character describing a sample in 'sample_column' for comparison
#' @param sample.2 Character describing a sample in 'sample_column' for comparison
#' @param minimum_count integer value specifying the minimum number of clones for a clonotype to be included in analysis
#' @param decreasing Boolean value specifying whether or not to sort by decreasing (or increasing) expansion
#' @param identity A value in 'clonotype_identity_mode' to specify clonotypes of interest (must first run the 'calculate_freq' function)
#' @param id.column A column name to search for identity. If NULL (default), will search 'clonotype_identity_mode'
#' @return a dataframe of clones sorted by their expansion.
#' @export
expanded_clones <- function(seurat_object,
                            sample.1=NULL,
                            sample.2=NULL,
                            minimum_count=0,
                            decreasing=TRUE,
                            identity=NULL,
                            id.column=NULL
                            ) {


  clonotype_column <- seurat_object@misc$clonotype_column
  sample_column <- seurat_object@misc$sample_column

  if (is.null(clonotype_column) | is.null(sample_column)) {
    stop("You must first perform 'calculate_freq'.")
  }
  
  # Determine column to search for identity
  if (is.null(id.column)==TRUE) {
    id.column <- "clonotype_identity_mode"
  } else {
    NULL # keep as provided value
  }

  ## Extract meta data
  # If calculate_freq performed with cell type provided
  if (is.null(seurat_object@misc$cell_type_column) == FALSE) { 
    clone_freqs <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,
                            c(clonotype_column, sample_column, "clonotype_count_per_sample","clonotype_freq_per_sample",id.column)]
  
  # If calculate_freq performed without cell type provided but alternative id.column provided
  } else if (is.null(id.column)==FALSE) {
    clone_freqs <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,
                                           c(clonotype_column, sample_column, "clonotype_count_per_sample","clonotype_freq_per_sample",id.column)]
    
  # If calculate_freq performed without cell type
  } else {
    clone_freqs <- seurat_object@meta.data[is.na(seurat_object@meta.data[clonotype_column])==FALSE,
                            c(clonotype_column, sample_column, "clonotype_count_per_sample","clonotype_freq_per_sample")]
  }

  # Get distinct rows (because of duplicates)
  clone_freqs <- distinct(clone_freqs)

  # Select clone identities of interest
  if (is.null(identity)==FALSE) {
    
    if (is.null(id.column)==FALSE) {
      
      clone_freqs <- subset(clone_freqs, clone_freqs[id.column]==identity)
    
    } else {
      
      clone_freqs <- subset(clone_freqs, clonotype_identity_mode==identity)
      
    }
   
  }


  ### No sample specified ###
  if (is.null(sample.1) & is.null(sample.2)) {

    # Sort dataframe
    output <- clone_freqs[with(clone_freqs, order(clonotype_freq_per_sample, decreasing=decreasing)),]

    # Exclude clones that are below minimum_count
    output <- output[output$clonotype_count_per_sample > minimum_count,]

    # Remove row names
    row.names(output) <- NULL

    return(output) # sorted dataframe of clonotype-sample combos


  ### Sample.1 specified ###
  } else if (is.null(sample.1)==FALSE) {

    # Select specified samples
    clone_freqs.1 <- clone_freqs[clone_freqs[,sample_column] %in% c(sample.1),]

    # Calculate freq as mean of multiple samples
    if (is.null(seurat_object@misc$cell_type_column) == FALSE) {

      clone_freqs.1 <- clone_freqs.1 %>%
        group_by(pick(c(clonotype_column,"clonotype_identity_mode"))) %>%
        summarise(clonotype_freq = sum(clonotype_freq_per_sample) / length(c(sample.1)), # mean calculated this way bc clones not in all samples
                  clonotype_count = sum(clonotype_count_per_sample))

    # Same as above but if identity mode not calculated
    } else {

      clone_freqs.1 <- clone_freqs.1 %>%
        group_by(pick(c(clonotype_column))) %>%
        summarise(clonotype_freq = sum(clonotype_freq_per_sample) / length(c(sample.1)), # mean calculated this way bc clones not in all samples
                  clonotype_count = sum(clonotype_count_per_sample))

    }

    ### Sample.2 is not specified ###
    if (is.null(sample.2)) {

      # Sort dataframe
      output <- clone_freqs.1[with(clone_freqs.1, order(clonotype_freq, decreasing=decreasing)),]

      # Exclude clones that are below minimum_count
      output <- output[output$clonotype_count > minimum_count,]

      # Remove row names
      row.names(output) <- NULL

      return(output) # sorted dataframe of clonotypes in specified sample (or aggregated samples)

    ### Sample.2 is specified ###
    } else {

      # Select specified samples
      clone_freqs.2 <- clone_freqs[clone_freqs[,sample_column] %in% c(sample.2),]

      # Calculate freq as mean of multiple samples
      clone_freqs.2 <- clone_freqs.2 %>%
        group_by(pick(c(clonotype_column))) %>%
        summarise(clonotype_freq = sum(clonotype_freq_per_sample) / length(c(sample.2)), # mean calculated this way bc clones not in all samples
                  clonotype_count = sum(clonotype_count_per_sample))


      # Merge clone freqs from samples 1 and 2
      output <- merge(clone_freqs.1, clone_freqs.2, by=c(clonotype_column),
                           all=TRUE, suffixes=c("_sample.1","_sample.2"))

      # Replace NA with 0
      output[is.na(output)] <- 0

      # Include total count
      output$combined_count <- output$clonotype_count_sample.1 + output$clonotype_count_sample.2

      # Calculate change in freq
      output$freq_change <- output$clonotype_freq_sample.1 - output$clonotype_freq_sample.2

      # Calculate fold change in freq
      output$freq_FC <- output$clonotype_freq_sample.1 / output$clonotype_freq_sample.2

      # Calculate change in count
      output$count_change <- output$clonotype_count_sample.1 - output$clonotype_count_sample.2

      # Calculate fold change in count
      output$count_FC <- output$clonotype_count_sample.1 / output$clonotype_count_sample.2

      # Sort dataframe
      output <- output[with(output, order(freq_FC, count_change, decreasing=decreasing)),]

      # Exclude clones that are below minimum_count
      output <- output[output$combined_count > minimum_count,]

      # Remove row names
      row.names(output) <- NULL

      return(output) # sorted dataframe of clonotypes expanded between sample.1 and sample.2

    }

  ### Only sample.2 is specified ###
  } else {

    stop("You cannot provide sample.2 without providing sample.1")

  }

}

