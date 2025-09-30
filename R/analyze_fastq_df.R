#' Analyze a FASTQ file and return a data.frame
#'
#' @param file_path Path to FASTQ file.
#' @param min_length Minimum read length (default: NULL = no filter).
#' @param min_avg_qual Minimum average base quality (default: NULL = no filter).
#' @param min_gc_content Minimum GC% (default: NULL = no filter).
#'
#' @return A data.frame with columns:
#'   - id
#'   - length
#'   - avg_quality
#'   - gc_content
#' @export
analyze_fastq_df <- function(file_path, min_length = NULL, min_avg_qual = NULL, min_gc_content = NULL) {
  res <- analyze_fastq_r(file_path, min_length, min_avg_qual, min_gc_content)
  as.data.frame(res, stringsAsFactors = FALSE)
}
