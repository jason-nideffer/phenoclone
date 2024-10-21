
#' Highlighting Multiple Clonotypes
#'
#' This function plots a UMAP where cells can be colored according to a variable.
#' Multiple clonotypes will also be highlighted.
#'
#' @param seurat_object seurat object to operate on
#' @param clonotypes Character vector specifying clonotypes within the 'clonotype_column'
#' @param highlight_colors Vector of colors to use for the highlighted clonotypes
#' @param other_color Color to use for the non-highlighted cells
#' @param highlight_size Number representing the point size of the highlighted clonotypes
#' @param other_size Number representing the point size of the non-highlighted cells
#' @param other_alpha Number (0-1) corresponding to the transparency of non-highlighted cells
#' @return a ggplot object that is plotted and can be stored/modified
#' @export
highlight_multiple_clonotypes <- function(seurat_object,
                                          clonotypes,
                                          highlight_colors = seurat_object@misc$colors,
                                          other_color = 'lightgray',
                                          highlight_size = 1,
                                          other_size = 0.3,
                                          other_alpha = 0.1) {

  # Retrieve clonotype column name
  clonotype_column <- seurat_object@misc$clonotype_column

  all_clones_data = data.frame()

  for (clonotype in clonotypes) {
    clone_data = subset(seurat_object@meta.data, seurat_object@meta.data[clonotype_column]==clonotype)

    clone_data = seurat_object@reductions$umap@cell.embeddings[rownames(clone_data),]
    clone_data = as.data.frame(clone_data)
    clone_data$clonotype <- clonotype

    all_clones_data <- rbind(clone_data, all_clones_data)
  }

  clone_plot <- geom_point(data=all_clones_data, aes(x=UMAP_1, y=UMAP_2,
                                                     group=clonotype, color=clonotype),
                           size=highlight_size, alpha=1)

  umap_data <- seurat_object@reductions$umap@cell.embeddings
  umap_data <- as.data.frame(umap_data)

  umap_plot <- geom_point(data=umap_data, aes(x=UMAP_1, y=UMAP_2),
                          size=other_size, alpha=other_alpha, color=other_color)

  plot <- ggplot() +
    umap_plot +
    clone_plot + scale_color_manual(values=highlight_colors, name=clonotype_column) +
    theme(
      # Hide panel borders and remove grid lines
      panel.border = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change axis line
      axis.line = element_line(colour = "black"),
      text = element_text(size = 15),
      # Remove gray background on legend
      legend.key=element_rect(fill="white")
    ) +
    guides(colour = guide_legend(override.aes = list(alpha=1, size=2)
    ))

  return(plot)
}
