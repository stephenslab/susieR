# Register downstream methods for the engine's exported and internal generics.
# No upstream namespace binding or existing method is replaced.
.engine <- function(name) get(name, envir = asNamespace("susieR"), inherits = FALSE)
.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  methods <- ls(ns, pattern = "\\.slide_individual$")
  for (method in methods) {
    generic <- sub("\\.slide_individual$", "", method)
    if (!exists(generic, envir = asNamespace("susieR"), inherits = FALSE))
      stop("This susieR installation lacks the required extension hook: ", generic)
    registerS3method(generic, "slide_individual", get(method, ns),
                     envir = asNamespace("susieR"))
  }
}
