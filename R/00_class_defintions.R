methods::setClass("CmdStanMCMC")
methods::setClass("CmdStanFit")
methods::setClass("R6")
methods::setClassUnion("stan_fit_classes", c("CmdStanMCMC", "CmdStanFit", "R6"))

methods::setClass(
  "mbsem",
  methods::representation(
    major_parameters = "data.frame",
    minor_factor_matrix = "data.frame",
    data_list = "list",
    stan_fit = "stan_fit_classes"
  )
)

new_mbsem <- function() {
  methods::new("mbsem")
}
