#' Interpret Functional Modules using LLM (Generic Dispatcher)
#'
#' @description
#' A generic dispatcher that routes to the appropriate LLM interpretation
#' function based on the class of `object`:
#' \itemize{
#'   \item `functional_module` class → \code{llm_interpret_functional_module()}
#'   \item `list` class (multi-omics object) → \code{llm_interpret_multi_omics_module()}
#' }
#' End users should always call \code{llm_interpret_module()} directly and
#' never need to call the underlying functions themselves.
#'
#' @param object Either a \code{functional_module} S4 object (single-omics) or
#'   a \code{list} object (multi-omics, produced by multi-omics pipeline).
#' @param ... Additional arguments passed to the underlying implementation
#'   function. See \code{\link{llm_interpret_functional_module}} or
#'   \code{\link{llm_interpret_multi_omics_module}} for the full parameter list.
#'
#' @return The updated object with LLM interpretation results stored in the
#'   appropriate slot (\code{llm_module_interpretation}).
#'
#' @examples
#' \dontrun{
#' # Single-omics: object is a functional_module class
#' result <- llm_interpret_module(
#'   object = my_functional_module,
#'   api_key = "your_api_key",
#'   embedding_output_dir = "/path/to/embedding/output"
#' )
#'
#' # Multi-omics: object is a list
#' result <- llm_interpret_module(
#'   object = multi_omics_modules,
#'   api_key = "your_api_key",
#'   embedding_output_dir = "/path/to/embedding/output",
#'   orgdb = org.Hs.eg.db
#' )
#' }
#'
#' @seealso
#' \code{\link{llm_interpret_functional_module}},
#' \code{\link{llm_interpret_multi_omics_module}}
#'
#' @export
llm_interpret_module <- function(object, ...) {
  UseMethod("llm_interpret_module")
}


#' @describeIn llm_interpret_module Method for \code{functional_module} S4 objects (single-omics).
#' @export
llm_interpret_module.functional_module <- function(object, ...) {
  llm_interpret_functional_module(object = object, ...)
}


#' @describeIn llm_interpret_module Method for \code{list} objects (multi-omics).
#' @export
llm_interpret_module.list <- function(object, ...) {
  llm_interpret_multi_omics_module(object = object, ...)
}


#' @describeIn llm_interpret_module Default method — raises an informative error
#'   for unsupported object types.
#' @export
llm_interpret_module.default <- function(object, ...) {
  stop(
    "No method defined for object of class '", paste(class(object), collapse = ", "), "'.\n",
    "  'llm_interpret_module()' supports:\n",
    "    - 'functional_module' (single-omics)\n",
    "    - 'list'              (multi-omics)\n",
    call. = FALSE
  )
}
