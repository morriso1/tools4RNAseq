#' Reading in csvs as list of character vectors.
#'
#' @param path_glob_str Glob to folder of csv files to be read.
#'
#' @return List of character vectors e.g. Gene sets
#' @export
#'
#' @examples
read_in_csvs_as_list_of_char <- function(path_glob_str) {
  gene_list = list()
  path_glob <- Sys.glob(path_glob_str)
  for (i in seq_along(path_glob)) {
    gene_list[[i]] <- readr::read_csv(path_glob[[i]]) %>%
      pull()
    names(gene_list)[i] <-
      stringr::str_split(path_glob[[i]], '/', simplify = TRUE) %>%
      last() %>% stringr::str_split("\\.", simplify = TRUE) %>% dplyr::first()
  }
  return(gene_list)
}
