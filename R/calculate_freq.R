
#' Calculate Clonotype Frequencies
#'
#' This function takes a seurat object with meta data that already includes
#' TCR annotations and calculates the frequencies of different clonotypes
#' on a per-sample basis.
#'
#' @param seurat_object seurat object to operate on
#' @param clonotype_column Character describing name of column to use for clone calling
#' @param sample_column Character describing name of column to use for sample discrimination
#' @param cell_type_column Character describing name of column to use for cell_type identification (could be and annotation or cluster)
#' @param exclude_na Boolean; if TRUE (default), cells without clonotype designation will be excluded from frequency calculations
#' @param percent Boolean; if TRUE, frequencies will be multiplied by 100 (default = FALSE)
#' @return seurat_object with meta data updated to include frequency data
#' @export
calculate_freq <- function(seurat_object,
                           clonotype_column,
                           sample_column,
                           cell_type_column=NULL,
                           subject_column=NULL,
                           exclude_na=TRUE,
                           percent=FALSE) {


  #### Wiping data in case function has already been performed ####
  seurat_object@meta.data <- seurat_object@meta.data %>%
    select(-any_of(c("clonotype.identity_count_per_sample",
                     "clonotype_count_per_sample",
                     "cells_per_sample",
                     "clonotype_freq_per_sample",
                     "clonotype.identity_freq_per_sample",
                     "clonotype_identity_mode"
                     )))

  seurat_object@misc$clonotype_column <- NULL
  seurat_object@misc$sample_column <- NULL
  seurat_object@misc$cell_type_column <- NULL
  seurat_object@misc$subject_column <- NULL
  ####



  # Store rownames and create new metadata object
  if (exclude_na) {
    meta <- seurat_object@meta.data[is.na(seurat_object@meta.data[[clonotype_column]])==FALSE,]
  } else {
    meta <- seurat_object@meta.data
  }

  # Calculate clonotype frequencies and total cells per sample
  to_add <- meta %>%
    group_by(pick(c(subject_column, sample_column, clonotype_column, cell_type_column))) %>%
    summarise(clonotype.identity_count_per_sample=n()) %>%
    group_by(pick(c(subject_column, sample_column, clonotype_column))) %>%
    mutate(clonotype_count_per_sample=sum(clonotype.identity_count_per_sample)) %>%
    group_by(pick(c(subject_column, sample_column))) %>%
    mutate(cells_per_sample=sum(clonotype.identity_count_per_sample))

  # Calculate clonotype_freq_per_sample
  to_add$clonotype_freq_per_sample <-  to_add$clonotype_count_per_sample / to_add$cells_per_sample
  
  if (percent) {
    to_add$clonotype_freq_per_sample <- to_add$clonotype_freq_per_sample * 100
  }

  # If identity column is provided
  if (is.null(cell_type_column)==FALSE) {

    to_add$clonotype.identity_freq_per_sample <-  to_add$clonotype.identity_count_per_sample / to_add$cells_per_sample
    
    if (percent) {
      to_add$clonotype.identity_freq_per_sample <- to_add$clonotype.identity_freq_per_sample * 100
    }

    # Calculate Clone Identity Modes (Identity calculated based on mean freq across all samples...)
    clonotype_identity_modes <- to_add %>%
      group_by(pick(c(subject_column, clonotype_column, cell_type_column))) %>%
      summarise(freq_clono_sum=sum(clonotype.identity_freq_per_sample)) %>%
      group_by(pick(subject_column, clonotype_column)) %>%
      slice_max(n=1, order_by = freq_clono_sum, na_rm = TRUE, with_ties = FALSE)

    clonotype_identity_modes <- select(clonotype_identity_modes, -c(freq_clono_sum))
    colnames(clonotype_identity_modes) <- c(subject_column, clonotype_column, "clonotype_identity_mode")

    to_add <- merge(to_add, clonotype_identity_modes, by=c(subject_column, clonotype_column))

  # If identity column is not provided
  } else {

    # Remove clonotype identity count if identity is not provided
    drops <- c("clonotype.identity_count_per_sample")
    to_add <- to_add[ , !(names(to_add) %in% drops)]

  }

  # Get rownames from original metadata
  old_meta <- seurat_object@meta.data
  rownames <- row.names(old_meta)

  # Assign order column to original metadata
  old_meta$order_column <- 1:nrow(old_meta)

  # Merge new data with original metadata
  new_meta <- merge(old_meta, to_add, by=c(subject_column, sample_column, clonotype_column, cell_type_column), all=TRUE)

  # Restore metadata order
  new_meta <- new_meta[order(new_meta$order_column), ]

  # Remove order column
  drops <- c("order_column","cells_per_sample")
  new_meta <- new_meta[ , !(names(new_meta) %in% drops)]

  # Add back rownames
  row.names(new_meta) <- rownames

  # Reassign metadata
  seurat_object@meta.data <- new_meta

  # Store column names used for clonotype, samples, and identity
  seurat_object@misc$clonotype_column <- clonotype_column
  seurat_object@misc$sample_column <- sample_column
  seurat_object@misc$cell_type_column <- cell_type_column
  seurat_object@misc$subject_column <- subject_column

  return (seurat_object)
}
