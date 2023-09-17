# Load libraries  ---------------------------------------------------------
library(pacman)
p_load(googledrive, glue, purrr)

#--------------------------------------------
## Set paths
#--------------------------------------------
paths <- list(
  inputPath = file.path("inputs"),
  outputPath = file.path("outputs")
)

#Download data from google drive 
drive_url <- 'https://drive.google.com/drive/folders/1B745U29C-8PpFCySPCu2IQzfXGTvmdAf'
folder_id <- drive_get(as_id(drive_url))
files <- drive_ls(folder_id)

# Check if the local directory exists, and create it if it doesn't
if (!file.exists(outputPath)) {
  dir.create(outputPath, recursive = TRUE)
  cat("Local directory created:", outputPath, "\n")
}

purrr::walk2(files$name, files$id, ~drive_download(file = drive_get(.y), 
                                   path = file.path(outputPath, .x)))


