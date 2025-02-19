## code to prepare `exp_design_example` dataset goes here

# Code to generate this data is available at:
# https://github.com/l-a-yates/cpam_manuscript/data

temp_file <- tempfile(fileext = ".rds")
download.file(
  "https://raw.githubusercontent.com/l-a-yates/cpam_manuscript/main/data/exp_design_example.rds",
  destfile = temp_file
)
exp_design_example <- readRDS(temp_file)

usethis::use_data(exp_design_example, overwrite = TRUE)
