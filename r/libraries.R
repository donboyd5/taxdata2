
library(Rcpp)

library(tidyverse)
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library(readr)
library(vroom)
library(readxl)
library(arrow)
library(fs)

library(scales)
library(hms) # for times
library(lubridate) # for date/times
library(vctrs)
library(skimr)

library(grDevices)
library(knitr)
library(kableExtra)
library(ggrepel)

# remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")
library(tabulizer) # get data from pdfs

# packages for matching ----
library(StatMatch)
library(RANN)
library(rlemon)
# library(optmatch) # not needed?
# library(FNN) # get.knnx(data, query, k, algorithm) -- not needed?

# my packages ----
# renv::install(here::here("renv/cellar/btools_0.9.5.tar.gz"))
library(btools)
# renv::install(here::here("renv/cellar/bggtools_0.0.9200.tar.gz"))
library(bggtools)

library(DT) # for datatable

library(zoo) # for rollapply

library(bea.R)
library(fredr)
library(httr)
library(jsonlite)