.onLoad <- function(libname, pkgname) {
  # packageStartupMessage("Welcome dude!")
}

.onUnload <- function (libpath) {
  # packageStartupMessage("Let's git out of here!")
  library.dynam.unload("LDtools", libpath)
}