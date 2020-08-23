#' Wrapper function for Limma TopTable function. Returns vector of genes ids.
#'
#' @param p_val_cutoff
#' @param LFC_cutoff
#' @param gene_id_of_choice Can be flybase id, symbol etc.
#' @param direction up, down or both
#' @param ... Passing the dots to TopTable function.
#'
#' @return
#' @export
#'
#' @examples
differential_expressed_genes_per_contrast_char <-
  function(p_val_cutoff,
           LFC_cutoff,
           gene_id_of_choice = "flybase_ids",
           direction = 'Both',
           ...) {
    top_dif_genes_df <- topTable(...) %>% as_tibble(rownames = "flybase_ids")
    top_dif_genes_df <- top_dif_genes_df %>%
      filter(adj.P.Val < p_val_cutoff, abs(logFC) > LFC_cutoff)

    if (direction == "Up") {
      top_dif_genes_df <- top_dif_genes_df %>% filter(logFC > 0)
    } else if (direction == "Down") {
      top_dif_genes_df <- top_dif_genes_df %>% filter(logFC < 0)
    }

    gene_vec <- top_dif_genes_df %>% pull(gene_id_of_choice) %>%
      as.character()

    gene_vec <- gene_vec[!is.na(gene_vec)]

    return(gene_vec)

  }

#' Wrapper function for Limma TopTable function. Returns dataframe.
#'
#' @param p_val_cutoff
#' @param LFC_cutoff
#' @param ... Passing the dots to TopTable function.
#'
#' @return Data Frame of differentially expressed genes from toptable.
#' @export
#'
#' @examples
differential_expressed_genes_per_contrast_df <-
  function(p_val_cutoff,
           LFC_cutoff,
           ...) {
    top_dif_genes_df <- topTable(...) %>% as_tibble(rownames = "flybase_ids")
    gene_list_df <- top_dif_genes_df %>%
      filter(P.Value < p_val_cutoff, abs(logFC) > LFC_cutoff)

    return(gene_list_df)
  }

#' Multiple contrast wrapper to identify differentially expressed genes
#'
#' @param constrast_matrix
#' @param target_function "differential_expressed_genes_per_contrast_char" or "differential_expressed_genes_per_contrast_df"
#' @param annotation_for_list_names
#' @param ... Pass the dots to the target_function. This then passes to the TopTable function.
#'
#' @return list of character vectors or list of dataframes.
#' @export
#'
#' @examples
multiple_constrast_top_genes <-
  function(constrast_matrix,
           target_function,
           annotation_for_list_names = '',
           ...) {
    multiple_gene_lists = list()
    i <- 1
    for (n_constrast in constrast_matrix %>% colnames()) {
      multiple_gene_lists[[i]] <- target_function(coef = n_constrast,
                                                  ...)
      names(multiple_gene_lists)[i] <- paste0(n_constrast, annotation_for_list_names)
      i <- i + 1
    }
    return(multiple_gene_lists)
  }
