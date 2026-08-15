## install_smcure.R
##
## Installs the `smcure` package (Fit Semiparametric Mixture Cure Models),
## the software the manuscript reports using (Section 2.5.3). smcure is on
## CRAN, but if install.packages("smcure") is unavailable in your network
## environment, this installs it directly from its official CRAN mirror on
## GitHub (https://github.com/cran/smcure), which only requires `survival`,
## `stats`, and `graphics` (all base/recommended R packages).

if (!requireNamespace("smcure", quietly = TRUE)) {
  ok <- tryCatch({
    install.packages("smcure")
    requireNamespace("smcure", quietly = TRUE)
  }, error = function(e) FALSE)

  if (!isTRUE(ok)) {
    message("install.packages('smcure') failed or CRAN unreachable; ",
            "falling back to the CRAN mirror on GitHub...")
    tmp <- tempfile(fileext = ".tar.gz")
    download.file("https://codeload.github.com/cran/smcure/tar.gz/refs/heads/master",
                   destfile = tmp, mode = "wb")
    tmpdir <- tempfile()
    dir.create(tmpdir)
    utils::untar(tmp, exdir = tmpdir)
    pkg_dir <- list.files(tmpdir, full.names = TRUE)[1]
    install.packages(pkg_dir, repos = NULL, type = "source")
  }
}

library(smcure)
cat("smcure", as.character(utils::packageVersion("smcure")), "is installed and loaded.\n")
