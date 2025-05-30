#' #' Launch MAPA Shiny Application
#' #'
#' #' This function launches the MAPA Shiny application included in the 'mapa' package.
#' #' It locates the Shiny app directory within the installed package and runs the app.
#' #'
#' #' @return None.
#' #' @export
#' #'
#' #' @note This function does not take any arguments.
#' #'       Ensure that the 'mapa' package is correctly installed with the Shiny app included.
#'
#' launch_mapa_shiny <-
#'   function() {
#'     if (!requireNamespace("shiny", quietly = TRUE)) {
#'       stop("The 'shiny' package is required but not installed. Please install 'shiny' to use this function.",
#'            call. = FALSE)
#'     }
#'
#'     appDir <- system.file("shinyapp", package = "mapa")
#'     if (appDir == "") {
#'       stop("Shiny app not found in the package", call. = FALSE)
#'     }
#'
#'     shiny::runApp(appDir)
#'   }
