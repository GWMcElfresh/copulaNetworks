#' Save a checkpoint RDS file
#'
#' @param object R object to save.
#' @param out_dir Output directory.
#' @param filename File name (default `checkpoint.rds`).
#' @return Invisibly returns the full file path.
#' @export
save_checkpoint <- function(object, out_dir, filename = "checkpoint.rds") {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(out_dir, filename)
  saveRDS(object, path)
  message("Saved checkpoint: ", path)
  invisible(path)
}

#' Load a checkpoint RDS file
#'
#' @param out_dir Directory containing the checkpoint.
#' @param filename File name (default `checkpoint.rds`).
#' @return The loaded object.
#' @export
load_checkpoint <- function(out_dir, filename = "checkpoint.rds") {
  path <- file.path(out_dir, filename)
  if (!file.exists(path)) {
    stop("Checkpoint not found: ", path, call. = FALSE)
  }
  readRDS(path)
}
