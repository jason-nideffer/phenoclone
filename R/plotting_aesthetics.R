#' Plotting Aesthetics
#'
#' This function takes a stores plotting aesthetics in a seurat object for simple, reproducible plotting.
#'
#' @param seurat_object seurat object to operate on
#' @param cell_type_order vector that includes cell types in desired order
#' @param colors vector that includes colors in corresponding order with 'cell_type_order'
#' @return a seurat object with new stored plotting aesthetics.
#' @export
plotting_aesthetics <- function(seurat_object,
                                cell_type_order,
                                colors) {

  seurat_object@misc$cell_type_order <- cell_type_order
  seurat_object@misc$colors <- colors

  return(seurat_object)

}
