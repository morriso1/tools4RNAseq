#' Saving list of dataframes to csvs.
#'
#' @param list
#' @param additional_annotation_string
#' @param directory
#' @param ... Passing the dots write.table function.
#'
#' @return
#' @export
#'
#' @examples
write_list_2_csv_directory <- function(list,
                                       additional_annotation_string = "",
                                       directory,
                                       ...) {
  if (!dir.exists(directory)) {
    dir.create(directory)
  }

  for (i in seq_along(list)) {
    file_name = paste0(
      getwd(),
      '/',
      directory,
      '/',
      names(list)[i],
      '_',
      additional_annotation_string,
      ".csv"
    )
    write.table(list[[i]],
                file = file_name,
                sep = ',', ...)
  }
}


#' Saving list of dataframes to html datatables.
#'
#' @param list
#' @param additional_annotation_string
#' @param directory
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
saving_html_datatables_from_list_of_tbls <- function(list,
           additional_annotation_string,
           directory, ...) {
    if (!dir.exists(directory)) {
      dir.create(directory)
    }
    for (i in seq_along(list)) {
      file_name = paste0(
        getwd(),
        '/',
        directory,
        '/',
        names(list)[i],
        '_',
        additional_annotation_string,
        ".html"
      )
      data_table = datatable(list[[i]], ...)
      htmlwidgets::saveWidget(data_table, file = file_name)
      names(list)[i] %>% print()
    }
  }
