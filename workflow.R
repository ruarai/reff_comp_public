
library(tidyverse)
library(processx)
source("R/functions.R")

round_id <- "2026-11-11"

repo_dir <- "../forecast-hub"

check_and_log_git_status(round_id)

# Switch to main branch and pull latest changes
message(run("git", c("fetch"), wd = repo_dir, echo = TRUE)$stdout)
message(run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)$stdout)
message(run("git", c("pull"), wd = repo_dir, echo = TRUE)$stdout)

# Run the forecasts
message(run("sbatch", c(str_c("--export=ALL,round_id=", round_id), "run_forecasts.slurm"))$stdout)

# Wait for the above to complete
message(run("squeue", c("-u", "ruarai"))$stdout)

combine_forecasts(round_id)


# run("git", c("branch", "-D", "uom-reff_comp"), wd = repo_dir, echo = TRUE)

# Switch to model branch and push new changes
run("git", c("pull", "origin", "uom-reff_comp"), wd = repo_dir, echo = TRUE)
run("git", c("reset", "--hard", "origin/uom-reff_comp"), wd = repo_dir, echo = TRUE)

run("git", c("switch", "uom-reff_comp"), wd = repo_dir, echo = TRUE)

fs::dir_copy("forecast_outputs/uom-reff_comp/", "../forecast-hub/model-output/uom-reff_comp/", overwrite = TRUE)
run("git", c("add", "."), wd = repo_dir, echo = TRUE)
run("git", c("commit", "-m", str_c("Add uom-reff_comp forecasts for ", round_id)), wd = repo_dir, echo = TRUE)

run("git", c("fetch", "origin"), wd = repo_dir, echo = TRUE)
run("git", c("rebase", "origin/main"), wd = repo_dir, echo = TRUE)
run("git", c("push", "-f", "origin", "uom-reff_comp"), wd = repo_dir, echo = TRUE)




# Return to main branch
run("git", c("switch", "main"), wd = repo_dir, echo = TRUE)
message(run("git", c("pull"), wd = repo_dir, echo = TRUE)$stdout)

run("git", c("reset", "--hard", "origin/main"), wd = repo_dir, echo = TRUE)

source("reporting.R")
