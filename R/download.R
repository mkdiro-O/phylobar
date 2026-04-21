#' Read an RDS file from Zenodo
#'
#' Downloads an RDS file from a Zenodo record using \code{zen4R::download_zenodo}
#' and reads it. If the download fails, a message is emitted and the vignette
#' exits.
#'
#' @param doi The DOI of the Zenodo record (e.g., "10.5281/zenodo.17477876").
#' @param file The name of the RDS file within the record.
#' @return An R object from the RDS file.
#' @importFrom knitr knit_exit
#' @export
#' @examples
#' dat <- read_zenodo_rds("10.5281/zenodo.17477876", "su-2020.rds")
read_zenodo_rds <- function(doi, file) {
    tryCatch(
        zen4R::download_zenodo(doi = doi, path = tempdir(),
                               files = file, quiet = TRUE),
        error = \(e) {
            message("Unable to download data: ", conditionMessage(e))
            knit_exit()
        }
    )
    readRDS(file.path(tempdir(), file))
}
