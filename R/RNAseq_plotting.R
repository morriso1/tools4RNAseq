#' Creates a MDS or PCA plot by passing GLimma::glMDSPlot output to ggplot2
#'
#' @param loc_exp_lcpm Log transformed count matrix.
#' @param top Top number of genes to be analysed e.g. top 500 genes
#' @param dimensions PC dimensions to be plotted.
#' @param gene.selection "pairwise" for MDS. "common" for PCA.
#' @param factor_one Column of DGElist$samples that contains first factor to be labelled by color e.g. treatment
#' @param factor_two Column of DGElist$samples that contains second factor to be labelled by shape e.g. RNAi
#' @param Title
#' @param geom_size
#' @param shape_legend_color
#'
#' @return plot
#' @export
#'
#' @examples
create_mds_or_pca_plot <- function(loc_exp_lcpm,
                                   top = 500,
                                   dimensions = c(1, 2),
                                   gene.selection = 'pairwise',
                                   factor_one = exp_DGElist$samples$treatment,
                                   factor_two = exp_DGElist$samples$RNAi,
                                   Title = "MDS plot using avg(logFC) as the distance",
                                   geom_size = 2.5,
                                   shape_legend_color = 'grey') {
  exp_lcpm_mds <-
    Glimma::glMDSPlot(
      loc_exp_lcpm,
      launch = FALSE,
      top = top,
      gene.selection = gene.selection
    )

  to_plot <- tibble(
    Dim1 = exp_lcpm_mds$points[, dimensions[1]],
    Dim2 = exp_lcpm_mds$points[, dimensions[2]],
    factor_one = factor_one,
    factor_two = factor_two,
  )

  mds_var_per <-
    round(exp_lcpm_mds$eig / sum(exp_lcpm_mds$eig) * 100, 1)

  if (gene.selection == 'pairwise') {
    x_axes_label <- 'MDS1 -'
    y_axes_label <- 'MDS2 -'
  }

  if (gene.selection == 'common') {
    x_axes_label <- 'PC1 -'
    y_axes_label <- 'PC2 -'
  }

  my_plot <- ggplot(to_plot, aes(Dim1, Dim2)) +
    geom_point(aes(shape = factor_two, color = factor_one),
               size = geom_size,
               stroke = 1) +
    guides(shape = guide_legend(override.aes = list(color = shape_legend_color))) +
    xlab(paste(x_axes_label, mds_var_per[1], "%", sep = "")) +
    ylab(paste(y_axes_label, mds_var_per[2], "%", sep = "")) +
    ggtitle(Title)
}

#' Create a heatmap using the heatmap.2 function from gplots
#'
#' @param loc_exp_lcpm Log transformed count matrix.
#' @param groups Column of DGElist$samples that contains group info.
#' @param sam_ids Column of DGElist$samples that contains sam_id info.
#' @param gene_set_charV Gene set as a vector. Identifies need to match rownames of loc_exp_lcpm e.g. flybase_ids
#' @param reg_express REGEX expression to identify columns of groups to include in heatmap
#' @param col_cell Vector of colors to label each of the sam_id columns in the heatmap
#' @param ... Passing the dots to heatmap.2 function.
#'
#' @return plot
#' @export
#'
#' @examples
create_heatmap_from_lcpm <-
  function(loc_exp_lcpm,
           groups,
           sam_ids,
           gene_set_charV,
           reg_express = 'LacZRNAi',
           col_cell = NA,
           ...) {
    gene_set_exp_lcpm <-
      loc_exp_lcpm[rownames(loc_exp_lcpm) %in% gene_set_charV, grep(reg_express, groups)]

    of_interest <-
      groups[grep(reg_express, groups)] %>% droplevels()
    sam_ids <- sam_ids[grep(reg_express, groups)] %>% droplevels()
    new_order <- sam_ids[base::order(of_interest)]
    gene_set_exp_lcpm <- gene_set_exp_lcpm[, as.vector(new_order)]

    display_order <-
      as.character(of_interest[base::order(of_interest)])


    if (!is.character(col_cell)) {
      print("Please provide col_cell parameter as character vector.")
      return(display_order)
    } else if (length(col_cell) != length(display_order)) {
      print(base::paste(
        "col_cell parameter should be length",
        length(display_order)
      ))
      return(display_order)
    } else {
      plt <- heatmap.2(gene_set_exp_lcpm, ColSideColors = col_cell,
                       ...)
      return(plt)

      return(display_order)
    }
  }

#' Creates heatmap using pheatmap function. Can have two rows of labels.
#'
#' @param sample_df df from DGElist$samples
#' @param x_lcpm Count matrix. Usually log transformed.
#' @param gene_set_char Gene set as a vector Identifies need to match rownames of loc_exp_lcpm e.g. flybase_ids
#' @param filt_1_quo Column name in sample_df. Quosure. Do not quote.
#' @param filt_1_vec Vector containing elements to select from filt_1 column of sample_df
#' @param filt_2_quo Column name in sample_df. Quosure. Do not quote.
#' @param filt_2_vec Vector containing elements to select from filt_2 column of sample_df
#' @param var_to_groupby Normally DGElist$samples$group. Normally a Factor
#' @param ... Passing the dots to pheatmap function.
#'
#' @return
#' @export
#'
#' @examples
pheatmap_wrapper_function <-
  function(sample_df,
           x_lcpm,
           gene_set_char,
           filt_1_quo =  RNAi,
           filt_1_vec = c('LacZRNAi'),
           filt_2_quo = treatment,
           filt_2_vec = NA,
           var_to_groupby = group,
           ...) {
    filt_1_quo <- rlang::enquo(filt_1_quo)
    filt_2_quo <- rlang::enquo(filt_2_quo)
    var_to_groupby <- rlang::enquo(var_to_groupby)

    desired_df <- sample_df %>%
      as_tibble(rownames = "SAM_ids") %>%
      arrange(!!var_to_groupby)

    if (all(is.na(filt_1_vec)) & all(is.na(filt_2_vec))) {
      print("must provide at least one filter.")
      return()
    }
    else if (all(!is.na(filt_1_vec)) & all(is.na(filt_2_vec))) {
      desired_df <- desired_df %>%
        filter(!!filt_1_quo %in% filt_1_vec) %>%
        dplyr::select(SAM_ids,!!filt_1_quo) %>%
        column_to_rownames(var = "SAM_ids")
    }
    else if (all(is.na(filt_1_vec)) & all(!is.na(filt_2_vec))) {
      desired_df <- desired_df %>%
        filter(!!filt_2_quo %in% filt_2_vec) %>%
        dplyr::select(SAM_ids,!!filt_2_quo) %>%
        column_to_rownames(var = "SAM_ids")
    }
    else {
      desired_df <- desired_df %>%
        filter(!!filt_1_quo %in% filt_1_vec, !!filt_2_quo %in% filt_2_vec) %>%
        dplyr::select(SAM_ids,!!filt_1_quo, !!filt_2_quo) %>%
        column_to_rownames(var = "SAM_ids")
    }

    plt <- x_lcpm %>% as_tibble(rownames = "Flybase_IDs") %>%
      dplyr::select(Flybase_IDs, rownames(desired_df)) %>%
      filter(Flybase_IDs %in% gene_set_char) %>%
      column_to_rownames(var = "Flybase_IDs") %>%
      pheatmap(
        scale = 'row',
        cluster_cols = FALSE,
        col = bluered(100),
        annotation_col = desired_df,
        border_color = NA,
        cex = 1,
        fontsize_row = 6,
        ...
      )
  }

#' Scale and tidy data prior to ggplot2 plot e.g. violin plot.
#'
#' @param loc_exp_lcpm Count matrix. Usually log transformed.
#' @param groups Column of DGElist$samples that contains group info.
#' @param RNAi Column of DGElist$samples that contains RNAi info.
#' @param treatment Column of DGElist$samples that contains treatment info.
#' @param sam_ids Column of DGElist$samples that contains sam_ids info.
#' @param gene_set_charV Gene set as a vector. Identifies need to match rownames of loc_exp_lcpm e.g. flybase_ids
#'
#' @return tidy df
#' @export
#'
#' @examples
scale_extract_gene_set_and_tidy_lcpm_data <- function(loc_exp_lcpm,
                                                      groups,
                                                      RNAi,
                                                      treatment,
                                                      sam_ids,
                                                      gene_set_charV) {
  scaled_lcpm <- loc_exp_lcpm %>% t() %>% scale() %>% t()

  names(groups) <- sam_ids
  names(RNAi) <- sam_ids
  names(treatment) <- sam_ids

  scaled_lcpm_gene_set <-
    scaled_lcpm[rownames(scaled_lcpm) %in% gene_set_charV,] %>%
    as_tibble(rownames = 'flybase_ids') %>%
    gather(SAM_ID, expression, -1)

  scaled_lcpm_gene_set <- scaled_lcpm_gene_set %>%
    mutate(
      group = recode(scaled_lcpm_gene_set$SAM_ID, !!!groups),
      treatment = recode(scaled_lcpm_gene_set$SAM_ID, !!!treatment),
      RNAi = recode(scaled_lcpm_gene_set$SAM_ID, !!!RNAi)
    )

  return(scaled_lcpm_gene_set)
}



#' Scale and tidy data prior to ggplot2 plot e.g. violin plot. Only using group info.
#'
#' @param loc_exp_lcpm Log transformed count matrix.
#' @param groups Column of DGElist$samples that contains group info.
#' @param sam_ids Column of DGElist$samples that contains sam_ids info.
#' @param gene_set_charV Gene set as a vector. Identifies need to match rownames of loc_exp_lcpm e.g. flybase_ids
#'
#' @return tidy_df
#' @export
#'
#' @examples
scale_extract_gene_set_and_tidy_lcpm_data_simple <-
  function(loc_exp_lcpm,
           groups,
           sam_ids,
           gene_set_charV) {
    scaled_lcpm <- loc_exp_lcpm %>% t() %>% scale() %>% t()

    names(groups) <- sam_ids

    scaled_lcpm_gene_set <-
      scaled_lcpm[rownames(scaled_lcpm) %in% gene_set_charV, ] %>%
      as_tibble(rownames = 'flybase_ids') %>%
      gather(SAM_ID, expression,-1)

    scaled_lcpm_gene_set <- scaled_lcpm_gene_set %>%
      mutate(group = recode(scaled_lcpm_gene_set$SAM_ID,!!!groups))

    return(scaled_lcpm_gene_set)
  }


#' Violin plot using tidy data
#'
#' @param tidy_gene_set_data
#'
#' @return
#' @export
#'
#' @examples
violinplot_from_tidy_gene_set_data <- function(tidy_gene_set_data) {
  ggplot(tidy_gene_set_data,
         aes(x = treatment, y = expression, fill = treatment)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    facet_grid(. ~ RNAi)
}


#' Save dotplots and csvs from reactome pathway analysis
#'
#' @param flybase_multiple_gene_list Gene list of flybase ids
#' @param entrez_ids Column of DGElist$genes that contains entrez id info.
#' @param start_of_file_name
#' @param save_directory
#' @param number_cat_dotplot
#' @param font_size_dotpot
#' @param height_dotplot
#' @param width_dotplot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
multiple_reactomePA_save_dotplot_and_csv <-
  function(flybase_multiple_gene_list,
           entrez_ids,
           start_of_file_name,
           save_directory,
           number_cat_dotplot = 15,
           font_size_dotpot = 8,
           height_dotplot = 10,
           width_dotplot = 15,
           ...) {
    if (!dir.exists(file.path(getwd(), save_directory, start_of_file_name))) {
      dir.create(file.path(getwd(), save_directory, start_of_file_name))
    }
    entrez_gl <- entrez_ids %>%  as_tibble(rownames = 'flybase_ids')

    for (i in seq_along(flybase_multiple_gene_list)) {
      gene_list <-
        tibble(flybase_ids = flybase_multiple_gene_list[[i]]) %>%
        left_join(entrez_gl) %>% drop_na() %>% pull(value) %>% as.character()
      enrichment_result <-
        enrichPathway(gene = gene_list, organism = 'fly', ...)
      complete_file_name <-
        file.path(
          getwd(),
          save_directory,
          start_of_file_name,
          paste0(
            start_of_file_name,
            '_',
            names(flybase_multiple_gene_list)[i]
          )
        )
      print(complete_file_name)
      enrichment_result@result %>% readr::write_csv(paste0(complete_file_name, '.csv'))
      dotplot(enrichment_result,
              showCategory = number_cat_dotplot,
              font.size = font_size_dotpot)
      ggsave(paste0(complete_file_name, '.pdf'),
             width = width_dotplot,
             height = height_dotplot)
    }
  }
