# Read an RDS file from URL

Downloads an online RDS file into a tempfile and reads it. If the
download fails, a message is emitted and the vignette exits.

## Usage

``` r
read_url_rds(url)
```

## Arguments

- url:

  The URL to download.

## Value

An R object from the RDS file.

## Examples

``` r
dat <- read_url_rds("https://zenodo.org/records/17477876/files/su-2020.rds?download=1")
```
