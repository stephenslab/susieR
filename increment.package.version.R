# A short R script to increment the version number in the package
# DESCRIPTION file.
out <- readLines("DESCRIPTION")
i   <- which(substr(out,1,8)  == "Version:")
if (length(i) > 0) {

  # Get the current version number.
  x <- trimws(unlist(strsplit(out[i],":",fixed = TRUE))[2])
  cat("Current package version:",x,"\n")
  x      <- unlist(strsplit(x,".",fixed = TRUE))
  n      <- length(x)
  x[n]   <- as.numeric(x[n]) + 1
  x      <- paste(x,collapse = ".")
  cat("New package version:    ",x,"\n")
  out[i] <- paste("Version:",x)
}
writeLines(out,"DESCRIPTION")
