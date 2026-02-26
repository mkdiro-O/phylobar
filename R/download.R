#' Read an RDS file from URL
#'
#' Downloads an online RDS file into a tempfile and reads it. If the download
#' fails, a message is emitted and the vignette exits.
#'
#' @param url The URL to download.
#' @return An R object from the RDS file.
#' @importFrom knitr knit_exit
#' @export
#' @examples
#' dat <- read_url_rds("https://zenodo.org/records/17477876/files/su-2020.rds?download=1")
read_url_rds <- function(url) {
    f <- tempfile()
    tryCatch(
        download.file(url, f),
        error = \(e) {
            message("Unable to download data: ", conditionMessage(e))
            knit_exit()
        }
    )
    readRDS(f)
}
