#' Plot Clonal Frequencies
#'
#' This function plots the frequencies of clonotypes across samples. Cell type
#' proportions within a given clonotype will also be shown.
#'
#' @param seurat_object Seurat object to operate on
#' @param clonotypes Vector of clonotypes to include for plotting
#' @param samples Vector of samples to consider for plotting
#' @param cell_type_order Vector of unique values from cell type column to specify order
#' @param cell_type_colors Vector of colors to use for coloring cell types (corresponds to 'cell_type_order')
#' @param plot.type Character; 'alluvial' or 'circularbar' to specify plot type
#' @return a ggplot object that is plotted and can be stored/modified.
#' @export
plot_freq <- function(seurat_object,
                      clonotypes,
                      samples,
                      cell_type_order = seurat_object@misc$cell_type_order,
                      cell_type_colors = seurat_object@misc$colors,
                      plot.type = "alluvial") {

  clonotype_column <- seurat_object@misc$clonotype_column
  sample_column <- seurat_object@misc$sample_column
  cell_type_column <- seurat_object@misc$cell_type_column

  meta <- seurat_object@meta.data[,c(clonotype_column,
                                     sample_column,
                                     cell_type_column,
                                     "clonotype.identity_freq_per_sample")]

  # Get unique combos of clonotype, sample, cell_type, and clonotype.identity_freq_per_sample
  meta <- meta %>%
    group_by(pick(c(clonotype_column,sample_column,cell_type_column,"clonotype.identity_freq_per_sample"))) %>%
    summarise()

  # Re-calculate clonotype.identity_freq_per_sample in case sample classification has changed
  meta <- meta %>%
    group_by(pick(c(clonotype_column,sample_column,cell_type_column))) %>%
    summarise(clonotype.identity_freq_per_sample = mean(clonotype.identity_freq_per_sample))

  # Subset on clonotypes
  meta <- subset(meta, meta[[clonotype_column]] %in% clonotypes)

  # Subset on samples
  meta <- subset(meta, meta[[sample_column]] %in% samples)

  # Make zero_template
  # Fill in zeros for clones not present in certain samples
  if (is.null(cell_type_order)==FALSE) {
    zero_template <- expand.grid(unique(meta[[clonotype_column]]),
                                 samples,
                                 cell_type_order)
  } else {
    zero_template <- expand.grid(unique(meta[[clonotype_column]]),
                                 samples,
                                 unique(meta[[cell_type_column]]))
  }

  colnames(zero_template) <- c(clonotype_column,sample_column,cell_type_column)

  # Add data on top of template
  meta <- merge(meta, zero_template, by=c(clonotype_column,sample_column,cell_type_column), all.y=TRUE)
  meta[is.na(meta)] <- 0

  # Add mode identity for sorting
  if (is.null(cell_type_order)==FALSE) {
    if("clonotype_identity_mode" %in% colnames(seurat_object@meta.data)) {

      mode_ident <- seurat_object@meta.data[,c(clonotype_column,"clonotype_identity_mode")]

      mode_ident <- mode_ident %>%
        group_by(pick(c(clonotype_column,"clonotype_identity_mode"))) %>%
        summarise() %>%
        group_by(pick(c(clonotype_column))) %>% 
        slice_sample(n = 1)

      meta <- merge(meta, mode_ident)

      # Tertiary sort based on cell type
      meta <- meta %>%
        arrange(factor(meta[,cell_type_column], levels = cell_type_order))

      # Secondary sort based on individual clone identity
      meta <- meta %>%
        arrange(factor(meta[,clonotype_column], levels = clonotypes))

      # Primary sort based on mode identity
      meta <- meta %>%
        arrange(factor(clonotype_identity_mode, levels = cell_type_order))

    } else {
      stop("To order by celltype, you must specify an identity column when running the function 'calculate freq'")
    }
  }

  # Plotting Alluvial
  if (plot.type == "alluvial") {

    meta[(meta)==0] <- 0.00000000000000000000000000000001

    # Add alluvium variable
    meta$alluvium <- unite(meta, "alluvium",
                           clonotype_column, cell_type_column, sep='_')$alluvium

    number_of_clones = length(unique(meta[[clonotype_column]]))
    number_of_clusters = length(cell_type_order)

    # Ordering data for alluvial
    meta <- meta %>%
      mutate(alluvium = factor(alluvium,
                               levels = unique(meta$alluvium)))

    # Specifying colors
    num_diff_colors_for_clones = 3
    repeats_for_clones = number_of_clones %/% num_diff_colors_for_clones
    remainder_for_clones = number_of_clones %% num_diff_colors_for_clones

    ramp <- c("#D4D4D4","#B4B4B4","#909090")

    clone_colors <- c(
      rep(ramp, repeats_for_clones),
      ramp[0:remainder_for_clones]
    )

    colors <- c(clone_colors, cell_type_colors)

    # Specify Order
    order <- c(unique(meta[[clonotype_column]]), cell_type_order)


    # Generating Palette
    palette = c()

    for (n in 1:length(order)) {
      palette[order[n]] <- colors[n]
    }
    
    plot <- ggplot(meta,
                   aes_string(x = sample_column, y="clonotype.identity_freq_per_sample",
                              stratum = "alluvium", alluvium = "alluvium",
                              fill = cell_type_column, label = clonotype_column,
                              group = clonotype_column)) +

      geom_flow(aes_string(fill = factor(meta[,clonotype_column])), stat = "alluvium",
                lode.guidance = "frontback", alpha=0.5) +

      scale_fill_manual(values = palette, breaks = cell_type_order) +

      geom_stratum(color=NA) +

      theme(legend.position = NULL, panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),panel.background = element_blank(),
            axis.line = element_line(colour = "black"), axis.text = element_text(size = 12),
            axis.title = element_text(size=15)) +

      xlab("") +
      ylab("Frequency")

    return(plot)

    # Plotting Circular Bar
  } else if (plot.type == "circularbar") {

    # Generate dataframe to make white circle
    circle_data <- data.frame(clonotype = unique(meta[,clonotype_column]),
                              val = rep(-(max(meta$clonotype.identity_freq_per_sample)*0.5),
                                        times = length(unique(meta[,clonotype_column]))
                              )
    )


    plot <- ggplot() +

      geom_col(
        data = meta,
        aes_string(
          x = factor(meta[,clonotype_column], level=unique(meta[[clonotype_column]])),
          y = "clonotype.identity_freq_per_sample",
          fill = factor(meta[,cell_type_column], level=unique(meta[[cell_type_column]]))
        ),
        position="stack",
        show.legend = TRUE,
        alpha = .9
      ) +

      # Make it circular
      coord_polar() +

      # Relabel Clones
      scale_x_discrete(labels=(1:length(clonotypes))) +

      # Make space for inner circle
      scale_y_continuous(limits = c(-(max(meta$clonotype.identity_freq_per_sample)*0.5),NA)) +

      # Define theme
      theme_minimal() +
      facet_wrap(sample_column, nrow=1) +

      # Fill inner circle
      geom_col(
        data = circle_data,
        aes(
          x = clonotype,
          y = val
        ),
        fill="white",
        alpha=1
      ) +

      xlab("") +
      ylab("Frequency") +

      # Set color
      if (is.null(cell_type_colors) == FALSE) {
        scale_fill_manual(values=cell_type_colors)
      }

    return(plot)

  }

}
