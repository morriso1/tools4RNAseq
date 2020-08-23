#' Add extra gene info to DGElist object.
#'
#' @param exp_DGElist DGElist object from edgeR package
#' @param org.Dm.eg.db
#' @param desired_columns
#' @param extra_gene_info
#' @param extra_gene_type_name
#' @param return_genes_info_only
#'
#' @return
#' @export
#'
#' @import org.Dm.eg.db
#'
#' @examples
edit_DGElist_genes <-
  function(exp_DGElist,
           org.Dm.eg.db,
           desired_columns = c('FLYBASE', 'SYMBOL', 'ONTOLOGY', 'PATH', 'ENTREZID'),
           extra_gene_info = NA,
           extra_gene_type_name = NA,
           return_genes_info_only = TRUE) {
    genes <- list()
    for (i in seq_along(desired_columns)) {
      genes[[i]] <- mapIds(
        org.Dm.eg.db,
        keys = rownames(exp_DGElist[["counts"]]),
        column = desired_columns[i],
        keytype = "ENSEMBL",
        multiVals = "first"
      )
      names(genes)[i] <- desired_columns[i]
    }
    genes <- as.data.frame(do.call(cbind, genes))

    if (!is.na(extra_gene_info[1])) {
      genes[[extra_gene_type_name]] <- extra_gene_info
    }
    if (return_genes_info_only){
      return(genes)
    } else{
      exp_DGElist$genes <- genes
      return(exp_DGElist)
    }

  }
