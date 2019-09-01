#' @title Startup Message
#' @description ARSeq- An automated RNASeq analysis pipeline
#' @param libname Library Name
#' @param pkgname Package name
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("#################################################################",
                              "Thank you for using ARSeq- An automated RNASeq analysis pipeline.",
                              "Check [[ https://ajitjohnson.com/arseq ]] for a complete tutorial",
                              "Hope you enjoy it- Ajit Johnson Nirmal (twitter: @ajitjohnson_n)",
                              "#################################################################",sep="\n"))
}
