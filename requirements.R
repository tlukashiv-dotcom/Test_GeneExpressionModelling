#' @title Project Dependencies Loader
#' @description Checks, optionally installs, and loads dependencies required
#' for the Embryo Molecule Counts modelling, validation, and visualization pipeline.

required_packages <- c(
  "mclust",
  "transport",
  "modeest",
  "mixtools",
  "mixdist",
  "DescTools",
  "ggplot2",
  "dplyr",
  "tidyr",
  "patchwork",
  "ggrepel",
  "foreach",
  "doParallel"
)

install_dependencies <- function(pkgs, install_missing = interactive()) {
  
  installed_pkgs <- rownames(installed.packages())
  missing_pkgs <- setdiff(pkgs, installed_pkgs)
  
  if (length(missing_pkgs) > 0) {
    
    if (isTRUE(install_missing)) {
      message("Installing missing dependencies: ", paste(missing_pkgs, collapse = ", "))
      install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
    } else {
      stop(
        "Missing required packages: ",
        paste(missing_pkgs, collapse = ", "),
        "
Install them manually or run install_dependencies(required_packages, install_missing = TRUE)."
      )
    }
  }
  
  invisible(lapply(pkgs, require, character.only = TRUE))
}

message("Checking and loading project dependencies...")
install_dependencies(required_packages)

message("--- Environment Summary ---")
message("R Version: ", R.version.string)
message("Platform : ", R.version$platform)
message("---------------------------")
