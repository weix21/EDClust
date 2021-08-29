#' Set up Julia
#'
#' This function initializes Julia and the EDClust.jl package.
#' The first time could be long due to precompilation.
#' Additionally, this function will call Julia and install the required packages
#' if they are missing.
#'
#' @param path path of Julia
#'
#' @return The Julia interface, which is an environment with the necessary methods provided by JuliaCall
#'
#' @export
#'
#' @import JuliaCall
#'
#' @examples
#'
#' julia <- julia_setup()
#'
#' julia <- julia_setup(path = "your Julia path")
#'
#'
setup_julia <- function(path = NA) {
  ## `RCall` needs to be precompiled with the current R.
  if (is.na(path)) {
    julia <- julia_setup(installJulia = TRUE)
  } else {
    julia <- julia_setup(JULIA_HOME = path)
  }
  julia$install_package_if_needed("Distributions")
  julia$install_package_if_needed("SpecialFunctions")
  julia$install_package_if_needed("StatsBase")
  if (julia$installed_package("EDClust") == "nothing") {
    julia_install_package("https://github.com/weix21/EDClust.jl.git")
  }
  julia$library("EDClust")
  return(julia)
}
