source("R/functions.R")
source("R/paper/common.R")

library(tinytable)


holiday_draws <- read_parquet("data/paper/holiday_draws.parquet")


plot_data <- holiday_draws |> 
  mutate(school_holiday_name = c("summer", "betweenterm")[ix_s]) |> 
  filter(location %in% locations_paper) |> 
  
  label_pathogen()


format_val <- function(x) {
  sprintf("%.2f",  exp(x))
}

table_data <- plot_data |> 
  group_by(school_holiday_name, location, pathogen) |> 
  summarise(
    median_Rt = median(school_holiday_Rt_effect_raw),
    lower_95_Rt = quantile(school_holiday_Rt_effect_raw,     0.05 / 2),
    upper_95_Rt = quantile(school_holiday_Rt_effect_raw, 1 - 0.05 / 2),
    
    Rt_CrI = glue::glue("{format_val(median_Rt)} ({format_val(lower_95_Rt)}–{format_val(upper_95_Rt)})"),
    
    bold_Rt = upper_95_Rt < 0,
    
    median_asc = median(school_holiday_asc_effect_raw),
    lower_95_asc = quantile(school_holiday_asc_effect_raw,     0.05 / 2),
    upper_95_asc = quantile(school_holiday_asc_effect_raw, 1 - 0.05 / 2),
    
    asc_CrI = glue::glue("{format_val(median_asc)} ({format_val(lower_95_asc)}–{format_val(upper_95_asc)})"),
    
    bold_asc = upper_95_asc < 0,
  ) |> 
  
  select(school_holiday_name, location, pathogen, Rt_CrI, asc_CrI, bold_Rt, bold_asc) |> 
  
  pivot_wider(
    names_from = school_holiday_name,
    values_from = c(Rt_CrI, asc_CrI, bold_Rt, bold_asc)
  ) |> 
  
  ungroup() |> 
  
  label_pathogen() |> 
  
  select(location, pathogen_label, Rt_CrI_betweenterm, asc_CrI_betweenterm,
         Rt_CrI_summer, asc_CrI_summer,
         bold_Rt_betweenterm, bold_Rt_summer,
         bold_asc_betweenterm, bold_asc_summer)

table_data |> 
  select(location, pathogen_label, Rt_CrI_betweenterm, asc_CrI_betweenterm,
         Rt_CrI_summer, asc_CrI_summer) |> 
  mutate(location = if_else(lag(location, default = "") != location &
                              lead(location, default = "") == location, location, "")) |> 
  
  `colnames<-`(c("Jurisdiction", "Pathogen", 
                 "$R(t)$ effect", "Reporting effect",
                 "$R(t)$ effect", "Reporting effect")) |> 
  
  tt() |> 
  
  group_tt(j = list("Between-term holidays" = 3:4, "Summer holidays" = 5:6)) |> 
  
  style_tt(j = 3, i = which(table_data$bold_Rt_betweenterm), bold = TRUE) |> 
  style_tt(j = 4, i = which(table_data$bold_asc_betweenterm), bold = TRUE) |> 
  style_tt(j = 5, i = which(table_data$bold_Rt_summer), bold = TRUE) |> 
  style_tt(j = 6, i = which(table_data$bold_asc_summer), bold = TRUE) |> 
  
  save_tt("paper_figures/supp/t_estimated_effects.tex", overwrite = TRUE)
