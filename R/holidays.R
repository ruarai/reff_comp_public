

get_school_holidays <- function() {
  # TODO verify dates here
  tribble(
    ~location, ~start_date,  ~end_date, ~holiday_name,
    
    "NSW", "2022-12-21", "2023-01-29", "Summer 2022/23",
    "NSW", "2023-04-07", "2023-04-23", "Autumn 2023",
    "NSW", "2023-07-01", "2023-07-16", "Winter 2023",
    "NSW", "2023-09-23", "2023-10-08", "Spring 2023",
    "NSW", "2023-12-20", "2024-01-29", "Summer 2023/24",
    "NSW", "2024-04-13", "2024-04-28", "Autumn 2024",
    "NSW", "2024-07-06", "2024-07-21", "Winter 2024",
    "NSW", "2024-09-28", "2024-10-13", "Spring 2024",
    "NSW", "2024-12-21", "2025-02-05", "Summer 2024/25", 
    "NSW", "2025-04-12", "2025-04-29", "Autumn 2025",    
    "NSW", "2025-07-05", "2025-07-21", "Winter 2025",    
    "NSW", "2025-09-27", "2025-10-13", "Spring 2025",    
    "NSW", "2025-12-20", "2026-02-01", "Summer 2025/26", 
    "NSW", "2026-04-03", "2026-04-21", "Autumn 2026",    
    "NSW", "2026-07-04", "2026-07-20", "Winter 2026",    
    "NSW", "2026-09-26", "2026-10-12", "Spring 2026",    
    "NSW", "2026-12-18", "2027-01-28", "Summer 2026/27", 
    
    "VIC", "2022-12-21", "2023-01-29", "Summer 2022/23",
    "VIC", "2023-04-07", "2023-04-23", "Autumn 2023",
    "VIC", "2023-06-24", "2023-07-09", "Winter 2023",
    "VIC", "2023-09-16", "2023-10-01", "Spring 2023",
    "VIC", "2023-12-21", "2024-01-29", "Summer 2023/24",
    "VIC", "2024-03-29", "2024-04-14", "Autumn 2024",
    "VIC", "2024-06-29", "2024-07-14", "Winter 2024",
    "VIC", "2024-09-21", "2024-10-06", "Spring 2024",
    "VIC", "2024-12-21", "2025-01-28", "Summer 2024/25",
    "VIC", "2025-04-05", "2025-04-21", "Autumn 2025",
    "VIC", "2025-07-05", "2025-07-20", "Winter 2025",
    "VIC", "2025-09-20", "2025-10-05", "Spring 2025",
    "VIC", "2025-12-20", "2026-01-27", "Summer 2025/26",
    "VIC", "2026-04-03", "2026-04-19", "Autumn 2026",
    "VIC", "2026-06-27", "2026-07-12", "Winter 2026",
    "VIC", "2026-09-19", "2026-10-04", "Spring 2026",
    "VIC", "2026-12-19", "2027-01-27", "Summer 2026/27",
    
    "QLD", "2022-12-10", "2023-01-22", "Summer 2022/23",
    "QLD", "2023-04-01", "2023-04-16", "Autumn 2023",
    "QLD", "2023-06-24", "2023-07-09", "Winter 2023",
    "QLD", "2023-09-16", "2023-10-02", "Spring 2023",
    "QLD", "2023-12-09", "2024-01-21", "Summer 2023/24",
    "QLD", "2024-03-29", "2024-04-14", "Autumn 2024",
    "QLD", "2024-06-22", "2024-07-07", "Winter 2024",
    "QLD", "2024-09-14", "2024-09-29", "Spring 2024",
    "QLD", "2024-12-14", "2025-01-27", "Summer 2024/25",
    "QLD", "2025-04-05", "2025-04-21", "Autumn 2025",
    "QLD", "2025-06-28", "2025-07-13", "Winter 2025",
    "QLD", "2025-09-20", "2025-10-06", "Spring 2025",
    "QLD", "2025-12-13", "2026-01-26", "Summer 2025/26",
    "QLD", "2026-04-03", "2026-04-19", "Autumn 2026",
    "QLD", "2026-06-27", "2026-07-12", "Winter 2026",
    "QLD", "2026-09-19", "2026-10-05", "Spring 2026",
    "QLD", "2026-12-12", "2027-01-26", "Summer 2026/27",
    
    "SA",  "2022-12-17", "2023-01-29", "Summer 2022/23",
    "SA",  "2023-04-15", "2023-04-30", "Autumn 2023",
    "SA",  "2023-07-08", "2023-07-23", "Winter 2023",
    "SA",  "2023-09-30", "2023-10-15", "Spring 2023",
    "SA",  "2023-12-16", "2024-01-28", "Summer 2023/24",
    "SA",  "2024-04-13", "2024-04-28", "Autumn 2024",
    "SA",  "2024-07-06", "2024-07-21", "Winter 2024",
    "SA",  "2024-09-28", "2024-10-13", "Spring 2024",
    "SA",  "2024-12-14", "2025-01-27", "Summer 2024/25",
    "SA",  "2025-04-12", "2025-04-27", "Autumn 2025",
    "SA",  "2025-07-05", "2025-07-20", "Winter 2025",
    "SA",  "2025-09-27", "2025-10-12", "Spring 2025",
    "SA",  "2025-12-13", "2026-01-26", "Summer 2025/26", 
    "SA",  "2026-04-11", "2026-04-26", "Autumn 2026",
    "SA",  "2026-07-04", "2026-07-19", "Winter 2026",
    "SA",  "2026-09-26", "2026-10-11", "Spring 2026",
    "SA",  "2026-12-12", "2027-01-31", "Summer 2026/27", 
    
    "WA",  "2022-12-16", "2023-01-31", "Summer 2022/23",
    "WA",  "2023-04-07", "2023-04-23", "Autumn 2023",
    "WA",  "2023-07-01", "2023-07-16", "Winter 2023",
    "WA",  "2023-09-23", "2023-10-08", "Spring 2023",
    "WA",  "2023-12-15", "2024-01-30", "Summer 2023/24",
    "WA",  "2024-03-29", "2024-04-14", "Autumn 2024",
    "WA",  "2024-06-29", "2024-07-14", "Winter 2024",
    "WA",  "2024-09-21", "2024-10-06", "Spring 2024",
    "WA",  "2024-12-13", "2025-02-04", "Summer 2024/25",
    "WA",  "2025-04-12", "2025-04-27", "Autumn 2025",
    "WA",  "2025-07-05", "2025-07-20", "Winter 2025",
    "WA",  "2025-09-27", "2025-10-12", "Spring 2025",
    "WA",  "2025-12-19", "2026-02-01", "Summer 2025/26",
    "WA",  "2026-04-03", "2026-04-19", "Autumn 2026",
    "WA",  "2026-07-04", "2026-07-19", "Winter 2026",
    "WA",  "2026-09-26", "2026-10-11", "Spring 2026",
    "WA",  "2026-12-18", "2027-01-31", "Summer 2026/27",
    
    "TAS", "2022-12-22", "2023-02-07", "Summer 2022/23",
    "TAS", "2023-04-06", "2023-04-25", "Autumn 2023",
    "TAS", "2023-07-08", "2023-07-24", "Winter 2023",
    "TAS", "2023-09-30", "2023-10-15", "Spring 2023",
    "TAS", "2023-12-22", "2024-02-06", "Summer 2023/24",
    "TAS", "2024-04-12", "2024-04-28", "Autumn 2024",
    "TAS", "2024-07-06", "2024-07-22", "Winter 2024",
    "TAS", "2024-09-28", "2024-10-13", "Spring 2024",
    "TAS", "2024-12-20", "2025-02-04", "Summer 2024/25",
    "TAS", "2025-04-12", "2025-04-27", "Autumn 2025",
    "TAS", "2025-07-05", "2025-07-21", "Winter 2025",
    "TAS", "2025-09-27", "2025-10-12", "Spring 2025",
    "TAS", "2025-12-19", "2026-02-03", "Summer 2025/26",
    "TAS", "2026-04-18", "2026-05-03", "Autumn 2026",
    "TAS", "2026-07-11", "2026-07-26", "Winter 2026",
    "TAS", "2026-10-03", "2026-10-18", "Spring 2026",
    "TAS", "2026-12-19", "2027-02-02", "Summer 2026/27",
    
    "NT",  "2022-12-17", "2023-01-30", "Summer 2022/23",
    "NT",  "2023-04-07", "2023-04-16", "Autumn 2023",
    "NT",  "2023-06-24", "2023-07-17", "Winter 2023",
    "NT",  "2023-09-23", "2023-10-08", "Spring 2023",
    "NT",  "2023-12-16", "2024-01-29", "Summer 2023/24",
    "NT",  "2024-04-06", "2024-04-14", "Autumn 2024",
    "NT",  "2024-06-22", "2024-07-15", "Winter 2024",
    "NT",  "2024-09-21", "2024-10-06", "Spring 2024",
    "NT",  "2024-12-14", "2025-01-28", "Summer 2024/25",
    "NT",  "2025-04-05", "2025-04-13", "Autumn 2025",
    "NT",  "2025-06-21", "2025-07-14", "Winter 2025",
    "NT",  "2025-09-20", "2025-10-06", "Spring 2025",
    "NT",  "2025-12-12", "2026-01-27", "Summer 2025/26",
    "NT",  "2026-04-03", "2026-04-12", "Autumn 2026",
    "NT",  "2026-06-20", "2026-07-13", "Winter 2026",
    "NT",  "2026-09-19", "2026-10-04", "Spring 2026",
    "NT",  "2026-12-11", "2027-01-26", "Summer 2026/27",
    
    "ACT", "2022-12-17", "2023-01-27", "Summer 2022/23",
    "ACT", "2023-04-07", "2023-04-23", "Autumn 2023",
    "ACT", "2023-07-01", "2023-07-16", "Winter 2023",
    "ACT", "2023-09-23", "2023-10-08", "Spring 2023",
    "ACT", "2023-12-16", "2024-01-29", "Summer 2023/24",
    "ACT", "2024-04-13", "2024-04-28", "Autumn 2024",
    "ACT", "2024-07-06", "2024-07-21", "Winter 2024",
    "ACT", "2024-09-28", "2024-10-13", "Spring 2024",
    "ACT", "2024-12-18", "2025-02-02", "Summer 2024/25",
    "ACT", "2025-04-11", "2025-04-28", "Autumn 2025",
    "ACT", "2025-07-04", "2025-07-21", "Winter 2025",
    "ACT", "2025-09-26", "2025-10-13", "Spring 2025",
    "ACT", "2025-12-19", "2026-01-29", "Summer 2025/26",
    "ACT", "2026-04-03", "2026-04-20", "Autumn 2026",
    "ACT", "2026-07-04", "2026-07-20", "Winter 2026",
    "ACT", "2026-09-26", "2026-10-12", "Spring 2026",
    "ACT", "2026-12-19", "2027-01-31", "Summer 2026/27"
  ) |> 
    mutate(
      start_date = ymd(start_date),
      end_date = ymd(end_date),
      school_holiday_type = case_when(
        str_detect(holiday_name, "Summer") ~ 1,
        str_detect(holiday_name, "Autumn") ~ 2,
        str_detect(holiday_name, "Winter") ~ 3,
        str_detect(holiday_name, "Spring") ~ 4,
      )
    )
}

get_school_holidays_long <- function() {
  get_school_holidays() |> 
    rowwise() |> 
    mutate(date = list(seq(start_date, end_date, by = "days"))) |> 
    unnest(date) |> 
    select(location, date, school_holiday_type)
}

get_public_holidays <- function() {
  locations_aus <- c("VIC", "NSW", "QLD", "NT", "TAS", "SA", "WA", "ACT")
  
  years_ahead <- 4
  
  aus_consistent <- tribble(
    ~date_origin, ~holiday_name,
    "2023-01-01", "New Year's Day",
    "2023-01-26", "Australia Day",
    "2023-04-25", "ANZAC Day",
    "2023-12-25", "Christmas Day",
    "2023-12-26", "Boxing Day",
  ) |> 
    rowwise() |> 
    mutate(date = list(ymd(date_origin) + years(0:years_ahead))) |> 
    unnest(date) |> 
    select(-date_origin) |> 
    expand_grid(location = locations_aus)
  
  nz_consistent <- tribble(
    ~date_origin, ~holiday_name,
    "2023-01-01", "New Year's Day",
    "2023-01-02", "Day after New Year's Day",
    "2023-02-06", "Waitangi Day",
    "2023-04-25", "ANZAC Day",
    "2023-12-25", "Christmas Day",
    "2023-12-26", "Boxing Day",
  ) |> 
    rowwise() |> 
    mutate(date = list(ymd(date_origin) + years(0:years_ahead))) |> 
    unnest(date) |> 
    select(-date_origin) |> 
    mutate(location = "NZ")
  
  easter_holidays <- tribble(
    ~date, ~holiday_name,
    "2023-04-07", "Easter",
    "2023-04-08", "Easter",
    "2023-04-09", "Easter",
    "2023-04-10", "Easter",
    
    "2024-03-29", "Easter",
    "2024-03-30", "Easter",
    "2024-03-31", "Easter",
    "2024-04-01", "Easter",
    
    "2025-04-18", "Easter",
    "2025-04-19", "Easter",
    "2025-04-20", "Easter",
    "2025-04-21", "Easter",
    
    "2026-04-03", "Easter",
    "2026-04-04", "Easter",
    "2026-04-05", "Easter",
    "2026-04-06", "Easter",
  ) |> 
    mutate(date = ymd(date)) |> 
    expand_grid(location = c("NZ", locations_aus))
  
  nz_variable_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "NZ", "2023-06-05", "King's Birthday",
    "NZ", "2023-07-14", "Matariki",
    "NZ", "2023-10-23", "Labour Day",
    
    "NZ", "2024-06-03", "King's Birthday",
    "NZ", "2024-06-28", "Matariki",
    "NZ", "2024-10-28", "Labour Day",
    
    "NZ", "2025-06-02", "King's Birthday",
    "NZ", "2025-06-20", "Matariki",
    "NZ", "2025-10-27", "Labour Day",
    
    "NZ", "2026-06-01", "King's Birthday",
    "NZ", "2026-07-10", "Matariki",
    "NZ", "2026-10-26", "Labour Day"
  ) |> 
    mutate(date = ymd(date))
  
  vic_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "VIC", "2023-03-13", "Labour Day",
    "VIC", "2023-06-12", "King's Birthday",
    "VIC", "2023-09-29", "Friday before the AFL Grand Final",
    "VIC", "2023-11-07", "Melbourne Cup",
    
    "VIC", "2024-03-11", "Labour Day",
    "VIC", "2024-06-10", "King's Birthday",
    "VIC", "2024-09-27", "Friday before the AFL Grand Final",
    "VIC", "2024-11-05", "Melbourne Cup",
    
    "VIC", "2025-03-10", "Labour Day",
    "VIC", "2025-06-09", "King's Birthday",
    "VIC", "2025-09-26", "Friday before the AFL Grand Final",
    "VIC", "2025-11-04", "Melbourne Cup",
    
    "VIC", "2026-03-09", "Labour Day",
    "VIC", "2026-06-08", "King's Birthday",
    "VIC", "2026-09-25", "Friday before the AFL Grand Final",
    "VIC", "2026-11-03", "Melbourne Cup"
  ) |> 
    mutate(date = ymd(date))
  
  nsw_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "NSW", "2023-06-12", "King's Birthday",
    "NSW", "2023-10-02", "Labour Day",
    
    "NSW", "2024-06-10", "King's Birthday",
    "NSW", "2024-10-07", "Labour Day",
    
    "NSW", "2025-06-09", "King's Birthday",
    "NSW", "2025-10-06", "Labour Day",
    
    "NSW", "2026-06-08", "King's Birthday",
    "NSW", "2026-10-05", "Labour Day"
  ) |> 
    mutate(date = ymd(date))
  
  act_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "ACT", "2023-03-13", "Canberra Day",
    "ACT", "2023-05-29", "Reconciliation Day",
    "ACT", "2023-06-12", "King's Birthday",
    "ACT", "2023-10-02", "Labour Day",
    
    "ACT", "2024-03-11", "Canberra Day",
    "ACT", "2024-05-27", "Reconciliation Day",
    "ACT", "2024-06-10", "King's Birthday",
    "ACT", "2024-10-07", "Labour Day",
    
    "ACT", "2025-03-10", "Canberra Day",
    "ACT", "2025-06-02", "Reconciliation Day",
    "ACT", "2025-06-09", "King's Birthday",
    "ACT", "2025-10-06", "Labour Day",
    
    "ACT", "2026-03-09", "Canberra Day",
    "ACT", "2026-06-01", "Reconciliation Day",
    "ACT", "2026-06-08", "King's Birthday",
    "ACT", "2026-10-05", "Labour Day"
  ) |> 
    mutate(date = ymd(date))
  
  qld_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "QLD", "2023-05-01", "Labour Day",
    "QLD", "2023-10-02", "King's Birthday",
    
    "QLD", "2024-05-06", "Labour Day",
    "QLD", "2024-10-07", "King's Birthday",
    
    "QLD", "2025-05-05", "Labour Day",
    "QLD", "2025-10-06", "King's Birthday",
    
    "QLD", "2026-05-04", "Labour Day",
    "QLD", "2026-10-05", "King's Birthday"
  ) |> 
    mutate(date = ymd(date))
  
  nt_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "NT", "2023-05-01", "May Day",
    "NT", "2023-06-12", "King's Birthday",
    "NT", "2023-08-07", "Picnic Day",
    
    "NT", "2024-05-06", "May Day",
    "NT", "2024-06-10", "King's Birthday",
    "NT", "2024-08-05", "Picnic Day",
    
    "NT", "2025-05-05", "May Day",
    "NT", "2025-06-09", "King's Birthday",
    "NT", "2025-08-04", "Picnic Day",
    
    "NT", "2026-05-04", "May Day",
    "NT", "2026-06-08", "King's Birthday",
    "NT", "2026-08-03", "Picnic Day"
  ) |> 
    mutate(date = ymd(date))
  
  sa_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "SA", "2023-03-13", "Adelaide Cup Day",
    "SA", "2023-06-12", "King's Birthday",
    "SA", "2023-10-02", "Labour Day",
    
    "SA", "2024-03-11", "Adelaide Cup Day",
    "SA", "2024-06-10", "King's Birthday",
    "SA", "2024-10-07", "Labour Day",
    
    "SA", "2025-03-10", "Adelaide Cup Day",
    "SA", "2025-06-09", "King's Birthday",
    "SA", "2025-10-06", "Labour Day",
    
    "SA", "2026-03-09", "Adelaide Cup Day",
    "SA", "2026-06-08", "King's Birthday",
    "SA", "2026-10-05", "Labour Day"
  ) |> 
    mutate(date = ymd(date))
  
  tas_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "TAS", "2023-03-13", "Eight Hours Day",
    "TAS", "2023-06-12", "King's Birthday",
    
    "TAS", "2024-03-11", "Eight Hours Day",
    "TAS", "2024-06-10", "King's Birthday",
    
    "TAS", "2025-03-10", "Eight Hours Day",
    "TAS", "2025-06-09", "King's Birthday",
    
    "TAS", "2026-03-09", "Eight Hours Day",
    "TAS", "2026-06-08", "King's Birthday"
  ) |> 
    mutate(date = ymd(date))
  
  wa_holidays <- tribble(
    ~location, ~date,          ~holiday_name,
    "WA", "2023-03-06", "Labour Day",
    "WA", "2023-06-05", "Western Australia Day",
    "WA", "2023-09-25", "King's Birthday",
    
    "WA", "2024-03-04", "Labour Day",
    "WA", "2024-06-03", "Western Australia Day",
    "WA", "2024-09-23", "King's Birthday",
    
    "WA", "2025-03-03", "Labour Day",
    "WA", "2025-06-02", "Western Australia Day",
    "WA", "2025-09-29", "King's Birthday",
    
    "WA", "2026-03-02", "Labour Day",
    "WA", "2026-06-01", "Western Australia Day",
    "WA", "2026-09-28", "King's Birthday"
  ) |> 
    mutate(date = ymd(date))
  
  bind_rows(
    aus_consistent,
    nz_consistent,
    nz_variable_holidays,
    easter_holidays,
    vic_holidays,
    nsw_holidays,
    act_holidays,
    qld_holidays,
    nt_holidays,
    sa_holidays,
    tas_holidays,
    wa_holidays
  ) |> 
    arrange(location, date)
}
