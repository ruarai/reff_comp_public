
library(tidyverse)
library(processx)

round_id <- "2025-05-08"


repo_dir <- "../forecast-hub"

# Switch to main branch and pull latest changes
run("git", c("fetch"), wd = repo_dir, echo = TRUE)
run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)
run("git", c("pull"), wd = repo_dir, echo = TRUE)

# TODO test this
run("sbatch", c(str_c("--export=ALL,round_id=", round_id), "run_forecasts.slurm"))

# Wait for the above to complete

combine_forecasts(round_id)

# Switch to model branch and push new changes
run("git", c("switch", "uom-reff-composure"), wd = repo_dir, echo = TRUE)

run("git", c("add", "."), wd = repo_dir, echo = TRUE)
run("git", c("commit", "-m", str_c("Added forecasts for ", round_id)), wd = repo_dir, echo = TRUE)
run("git", c("push", "-u", "origin", "uom-reff-composure"), wd = repo_dir, echo = TRUE)


# Return to main branch
run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)