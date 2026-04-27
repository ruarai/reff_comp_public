
library(httr)

# Define Melbourne coordinates
lat <- -37.8136
lon <- 144.9631

# Calculate date range (last 5 years from today)
end_date <- Sys.Date()
start_date <- end_date - as.difftime(5 * 365.25, units = "days")

# Construct the API URL for the Open-Meteo Archive
# Variables: relative_humidity_2m_mean
url <- paste0(
  "https://archive-api.open-meteo.com/v1/archive?",
  "latitude=", lat,
  "&longitude=", lon,
  "&start_date=", start_date,
  "&end_date=", end_date,
  "&daily=relative_humidity_2m_mean",
  "&timezone=Australia%2FSydney"
)

# Execute the GET request
response <- GET(url)
data_json <- jsonlite::fromJSON(content(response, "text", encoding = "UTF-8"))

humidity_data <- as.data.frame(data_json$daily) %>%
  rename(
    date = time,
    avg_relative_humidity = relative_humidity_2m_mean
  ) %>%
  mutate(date = as.Date(date)) |> 
  filter(date >= ymd("2023-01-01")) |> 
  mutate(rel_hum_smooth = zoo::rollmean(avg_relative_humidity, k = 28, fill = NA))

model_data <- x |> filter(location == "VIC") |> 
  filter(day_of_year < 366) |> 
  expand_grid(date = ymd(c("2023-01-01"), ymd("2024-01-01"), ymd("2025-01-01"), ymd("2026-01-01"))) |> 
  mutate(date = date + day_of_year)

ggplot(humidity_data) +
  geom_line(aes(x = date, y = rel_hum_smooth)) +
  
  geom_line(aes(x = date, y = 10 * humidity_scaled + 70),
            model_data)





