#' Save a checkpoint RDS file
#'
#' @param object R object to save.
#' @param outDir Output directory.
#' @param filename File name (default `checkpoint.rds`).
#' @return Invisibly returns the full file path.
#' @export
SaveCheckpoint <- function(object, outDir, filename = "checkpoint.rds") {
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
  }
  file_path <- file.path(outDir, filename)
  saveRDS(object, file_path)
  message("Saved checkpoint: ", file_path)
  invisible(file_path)
}

#' Load a checkpoint RDS file
#'
#' @param outDir Directory containing the checkpoint.
#' @param filename File name (default `checkpoint.rds`).
#' @return The loaded object.
#' @export
LoadCheckpoint <- function(outDir, filename = "checkpoint.rds") {
  file_path <- file.path(outDir, filename)
  if (!file.exists(file_path)) {
    stop("Checkpoint not found: ", file_path, call. = FALSE)
  }
  readRDS(file_path)
}
