
library(ggtext)

plot_theme_paper <- list(
  theme_minimal(),
  theme(panel.grid = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        text = element_text(family = "Helvetica", colour = "black"),
        line = element_line(linewidth = 0.7),
        axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm"), colour = "black"),
        axis.text.y = element_text(margin = margin(r = 0.3, unit = "cm"), colour = "black"),
        plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"),
        axis.ticks.length=unit(-0.1, "cm"),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        strip.text.x = element_text(size = 11, hjust = 0),
        plot.title = element_markdown(),
        plot.subtitle = element_markdown()
  )
)

plot_theme_annotation <- theme(plot.title = element_markdown(family = "Helvetica"),
                               plot.margin = margin(t = 0.2, l = 0.2, b = 0.2, r = 0.2, unit = "cm"))

element_gridline <- element_line(colour = "grey70", linetype = "28", linewidth = 0.3)
element_facet_background <- element_rect(fill = "white", colour = "grey50", linewidth = 0.5)


label_pathogen <- function(tbl) {
  tbl |> 
    mutate(pathogen_label = factor(
      pathogen,
      c("RSV", "flu", "SARSCOV2"),
      labels = c("RSV" = "RSV", "flu" = "Influenza", "SARSCOV2" = "SARS-CoV-2")
    ))
}

locations_paper <- c("NSW", "QLD", "VIC", "NT", "ACT")

process_perf_data <- function(perf_data) {
  
  processed_perf <- perf_data |> 
    mutate(horizon = as.numeric(date - origin_date)) |> 
    
    # Drop dates with known issues
    anti_join(
      tribble(
        ~location, ~pathogen, ~date,
        "WA", "flu", ymd("2025-08-16"),
        "WA", "flu", ymd("2025-08-23"),
        "WA", "flu", ymd("2025-08-30"),
        "WA", "RSV", ymd("2025-08-16"),
        "WA", "RSV", ymd("2025-08-23"),
        "WA", "RSV", ymd("2025-08-30")
      ),
      by = join_by(location, pathogen, date)
    ) |> 
    
    filter(date >= ymd("2025-01-01"),
           ymd(round_id) < ymd("2025-11-27"),
           location %in% locations_paper)
  
  return(processed_perf)
}

read_perf <- function(output_dir, model_name = output_dir) {
  summ_perf_file <- str_c("saved_results/", output_dir, "_performance.parquet")
  
  if(fs::file_exists(summ_perf_file)) {
    perf_data <- read_parquet(summ_perf_file) |> 
      collect() |> 
      mutate(.model = model_name)
  } else {
    perf_data <- open_dataset(
      fs::dir_ls(str_c("saved_results/", output_dir), recurse = TRUE, glob = "*performance.parquet")
    ) |> 
      collect() |> 
      mutate(.model = model_name)
  }
  
  return(perf_data)
}
