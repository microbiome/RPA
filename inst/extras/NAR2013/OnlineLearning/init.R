  warning("experimental")

  if (batch.size < 3) {
    warning("Minimum batch.size is 3. Setting batch.size = 3.")
    batch.size <- 3
  }

  if (is.null(cel.files) && !is.null(cel.path)) {
    cel.files <- list.celfiles(cel.path, full.names = TRUE)
  }

  #################################################################
  
  # NOTE: list CEL files in random order to avoid biases!
  # FIXME: add this
  cel.files.sampled <- cel.files # sample(cel.files)

