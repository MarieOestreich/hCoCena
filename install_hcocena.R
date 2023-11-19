install.packages('remotes')
remotes::install_github("tpq/propr")
remotes::install_url(url = "https://cran.r-project.org/src/contrib/Archive/MCDA/MCDA_0.0.24.tar.gz")
remotes::install_url(url="https://github.com/MarieOestreich/hCoCena/archive/refs/heads/version_1.2.zip", subdir="hCoCena-r-package", dependencies=TRUE, upgrade="never")
