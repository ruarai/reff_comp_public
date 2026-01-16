library(bayesplot)
library(tidybayes)
library(tidyverse)
library(patchwork)


source("R/functions.R")

jobs |> 
  pmap(function(location, pathogen) {
    make_fit_plots("2025-10-09", pathogen, location)
  })






