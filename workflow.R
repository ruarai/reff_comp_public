
library(tidyverse)
library(processx)
library(gert)
source("R/functions.R")

round_id <- "2025-05-08"

is_clean <- nrow(git_status()) == 0
if (!is_clean) {
  print("Working directory is not clean. Don't continue until all changes have been committed.")
}


repo_dir <- "../forecast-hub"

# Switch to main branch and pull latest changes
run("git", c("fetch"), wd = repo_dir, echo = TRUE);
run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)
run("git", c("pull"), wd = repo_dir, echo = TRUE)

# Run the forecasts
run("sbatch", c(str_c("--export=ALL,round_id=", round_id), "run_forecasts.slurm"))

# Wait for the above to complete

combine_forecasts(round_id)


# Switch to model branch and push new changes
run("git", c("pull", "origin", "uom-reff_comp"), wd = repo_dir, echo = TRUE)
# run("git", c("reset", "--hard", "origin/uom-reff_comp"), wd = repo_dir, echo = TRUE)

run("git", c("switch", "uom-reff_comp"), wd = repo_dir, echo = TRUE)

fs::dir_copy("forecast_outputs/uom-reff_comp/", "../forecast-hub/model-output/uom-reff_comp/", overwrite = TRUE)
run("git", c("add", "."), wd = repo_dir, echo = TRUE)
run("git", c("commit", "-m", str_c("Added forecasts for ", round_id)), wd = repo_dir, echo = TRUE)
run("git", c("push", "-u", "origin", "uom-reff_comp"), wd = repo_dir, echo = TRUE)


# Return to main branch
run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)
