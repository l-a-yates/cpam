## code to prepare `cpo_example` dataset goes here

# Code to generate this data is available at:
# https://github.com/l-a-yates/cpam_manuscript/data

temp_file <- tempfile(fileext = ".rds")
download.file(
  "https://raw.githubusercontent.com/l-a-yates/cpam_manuscript/main/data/cpo_example.rds",
  destfile = temp_file
)
cpo_example <- readRDS(temp_file)

usethis::use_data(cpo_example, overwrite = TRUE)
