library(bayesplot)
library(tidybayes)
library(tidyverse)
library(patchwork)


source("R/functions.R")

jobs |> 
  pmap(function(location, pathogen) {
    make_fit_plots("2025-06-19", pathogen, location)
  })






