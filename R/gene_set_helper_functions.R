#' Multiple contrast roast wrapper function - not mroast
#'
#' @param v voomed DGElist
#' @param gene_set character vector or list of character vectors containing gene sets. Ids must match v e.g. Flybase Ids.
#' @param design design matrix
#' @param contr_matrix contrast matrix
#' @param nrot number of rotations to pass to roast function.
#'
#' @return df
#' @export
#'
#' @examples
multiple_constrast_roast <-
  function(v, gene_set, design, contr_matrix, nrot = 9999) {
    roast_list = list()
    for (i in 1:ncol(contr_matrix)) {
      roast_list[[i]] <- roast(
        v,
        gene_set,
        design = design,
        contrast = contr_matrix[, i],
        nrot = 9999
      )$p.value %>% as_tibble(rownames = "Direction")
      names(roast_list)[i] <- colnames(contr_matrix)[i]
      print(colnames(contr_matrix)[i])
    }
    roast_df <- bind_rows(roast_list, .id = 'constrast')
    return(roast_df)
  }



#' Extract gene set character vectors from list of character vectors.
#'
#' @param list_of_tbls e.g. from dplyr group_split function.
#' @param gene_id
#' @param gene_set_name
#'
#' @return list of gene set character vectors.
#' @export
#'
#' @examples
creating_gene_sets <-
  function(list_of_tbls,
           gene_id = "entrez_gene",
           gene_set_name = "gs_name") {
    list_gene_sets <- list()
    for (i in seq_along(list_of_tbls)) {
      list_gene_sets[[i]] <- list_of_tbls[[i]] %>%
        pull(gene_id) %>%
        as.character()

      names(list_gene_sets)[i] <-
        list_of_tbls[[i]][[1, gene_set_name]]
    }
    return(list_gene_sets)
  }



#' Multiple contrast wrapper for limma camera function.
#'
#' @param contr_matrix
#' @param ... Passing dots to limma camera function.
#'
#' @return List of camera tibbles. Can later bind in one tibble using dplyr::bind_rows(.id = "contrast").
#' @export
#'
#' @examples
multiple_constrast_camera <-
  function(contr_matrix, ...) {
    camera_list = list()

    for (i in 1:ncol(contr_matrix)) {
      camera_list[[i]] <-
        camera(contrast = contr_matrix[, i], ...) %>%
        as_tibble(rownames = "Pathway or GO term")
      names(camera_list)[i] <- colnames(contr_matrix)[i]
      print(colnames(contr_matrix)[i])
    }
    return(camera_list)
  }

#' Multiple contrast wrapper for mroast camera function.
#'
#' @param contr_matrix
#' @param ... Passing dots to limma mroast function.
#'
#' @return List of mroast tibbles. Can later bind in one tibble using dplyr::bind_rows(.id = "contrast")
#' @export
#'
#' @examples
multiple_constrast_mroast <-
  function(contr_matrix, ...) {
    roast_list = list()
    for (i in 1:ncol(contr_matrix)) {
      roast_list[[i]] <- mroast(contrast = contr_matrix[, i], ...) %>%
        as_tibble(rownames = 'gene_set')
      names(roast_list)[i] <- colnames(contr_matrix)[i]
      print(colnames(contr_matrix)[i])
    }
    return(roast_list)
  }

#' Convert character vector of flybase_ids to character vector gene symbols.
#'
#' @param flybase_id_char
#' @param gene_symbols_with_rownames
#'
#' @return
#' @export
#'
#' @examples
convert_flybase_id_to_gene_symbol <-
  function(flybase_id_char,
           gene_symbols_with_rownames) {
    left_join(
      tibble(flybase_id = flybase_id_char),
      as_tibble(gene_symbols_with_rownames, rownames = 'flybase_id'),
      by = "flybase_id"
    ) %>% pull(value) %>% as.character() %>% return()

  }
