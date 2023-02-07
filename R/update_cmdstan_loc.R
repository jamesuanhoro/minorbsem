#' Update CmdStan location
#'
#' @param loc (string) Folder CmdStan lives on system
#' @examples
#' update_cmdstan_loc(loc = "~/cmdstan/")
#' @export
update_cmdstan_loc <- function(loc = "") {
  writeLines(
    text = loc,
    con = system.file("cmdstan_loc", package = "minorbsem")
  )
}
