





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






