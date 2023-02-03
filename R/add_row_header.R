add_row_header <- function(kbl_object, table_to_print, search_term) {
  if (any(table_to_print[1] == search_term)) {
    kbl_object <- kbl_object |>
      kableExtra::pack_rows(
        search_term,
        which(table_to_print[1] == search_term)[1],
        rev(which(table_to_print[1] == search_term))[1]
      )
  }
  return(kbl_object)
}
