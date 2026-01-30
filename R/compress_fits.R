




compress_fit <- function(fit_file, out_file) {
  draws <- tidy_draws(read_rds(fit_file))
  
  target_schema <- schema(lapply(draws_reduced, function(x) {
    if (is.numeric(x) && !is.integer(x)) float32() else infer_type(Array$create(x))
  }))
  
  draws_reduced |>
    filter(.draw %% 4 == 0) |>
    as_arrow_table(schema = target_schema) |>
    write_parquet(out_file, compression = "zstd")
}

files <- fs::dir_ls("results", recurse = TRUE, glob = "*.rds")

for(fit_file in files) {
  print(fit_file)
  new_name <- str_c(
    fs::path_ext_remove(fit_file),
    "_reduced_draws.parquet"
  )
  
  compress_output(fit_file, new_name)
  fs::file_delete(fit_file)
}






