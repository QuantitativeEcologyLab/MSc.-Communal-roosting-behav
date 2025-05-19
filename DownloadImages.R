library(httr)
library(jsonlite)
library(tidyverse)

order_name <- "Strigiformes"
output_folder <- "phylopic_species_images"
dir.create(output_folder, showWarnings = FALSE)

trait_data <- read_csv("Chapters/Bird_data_clean.csv")

species_list <- trait_data %>%
  filter(Order == order_name) %>%
  pull(Species) %>%
  unique()

cat("Attempting to fetch images for", length(species_list), "species.\n")

download_phylopic_png <- function(species_name, outdir, verbose=FALSE) {
  file_name <- gsub(" ", "_", species_name)
  png_path <- file.path(outdir, paste0(file_name, ".png"))
  
  if (file.exists(png_path)) {
    if (verbose) cat(species_name, ": Already downloaded.\n")
    return("already_downloaded")
  }
  
  # Step 1: Search for matching images
  query_url <- paste0(
    "https://api.phylopic.org/v2/images?names=",
    URLencode(species_name),
    "&options=pngFiles"
  )
  if (verbose) cat("Query URL for", species_name, ":", query_url, "\n")
  
  resp <- try(httr::GET(query_url), silent=TRUE)
  if (inherits(resp, "try-error")) {
    if (verbose) cat("httr::GET error for", species_name, "\n")
    return("GET_error")
  }
  if (resp$status_code != 200) {
    if (verbose) cat("Non-200 status for", species_name, ":", resp$status_code, "\n")
    return(paste0("HTTP_", resp$status_code))
  }
  
  raw <- rawToChar(resp$content)
  if (verbose) cat("Raw API content (first 200 chars):", substr(raw, 1, 200), "\n")
  result <- try(fromJSON(raw), silent=TRUE)
  if (inherits(result, "try-error")) {
    if (verbose) cat("JSON parse error for", species_name, "\n")
    return("JSON_error")
  }
  if (length(result$`_links`$images) == 0) {
    if (verbose) cat("No images found for", species_name, "\n")
    return("not_found")
  }
  
  image_uuid <- sub(".*/", "", result$`_links`$images[[1]]$href)
  image_info_url <- paste0("https://api.phylopic.org/v2/images/", image_uuid)
  if (verbose) cat("Image info URL:", image_info_url, "\n")
  
  info_resp <- try(httr::GET(image_info_url), silent=TRUE)
  if (inherits(info_resp, "try-error")) {
    if (verbose) cat("httr::GET error (info) for", species_name, "\n")
    return("info_GET_error")
  }
  if (info_resp$status_code != 200) {
    if (verbose) cat("Non-200 status for image info for", species_name, ":", info_resp$status_code, "\n")
    return(paste0("info_HTTP_", info_resp$status_code))
  }
  
  img_info_raw <- rawToChar(info_resp$content)
  img_info <- try(fromJSON(img_info_raw), silent=TRUE)
  if (inherits(img_info, "try-error")) {
    if (verbose) cat("Image info JSON parse error for", species_name, "\n")
    return("info_JSON_error")
  }
  if (is.null(img_info$pngFiles) || length(img_info$pngFiles) == 0) {
    if (verbose) cat("No PNG files for", species_name, "\n")
    return("no_png")
  }
  
  png_url <- img_info$pngFiles[[length(img_info$pngFiles)]]$url
  if (verbose) cat("Downloading PNG from:", png_url, "\n")
  tryCatch({
    download.file(png_url, png_path, mode="wb", quiet=!verbose)
    Sys.sleep(0.5)
    if (verbose) cat("Downloaded:", png_path, "\n")
    return("downloaded")
  }, error = function(e) {
    if (verbose) cat("Download error for", species_name, "\n")
    return("download_failed")
  })
}

# Load rphylopic
library(rphylopic)
# Get a single image UUID for a species
uuid <- get_uuid(name = "Acanthorhynchus_tenuirostris")
# Get the image for that UUID
img <- get_phylopic(uuid = uuid)
# But multiple silhouettes can exist per species...
uuid <- get_uuid(name = "Canis lupus", n = 5)

# How do I pick?!
# It's difficult without seeing the image itself, let's use:
img <- pick_phylopic(name = "Acanthorhynchus_tenuirostris", n = 4, view = 4)

# Try the first 3 species with verbose output
cat("---- DEBUG: First 3 species ----\n")
for (i in 1:min(3, length(species_list))) {
  res <- download_phylopic_png(species_list[i], output_folder, verbose=TRUE)
  cat(species_list[i], "->", res, "\n\n")
}

# Then run the rest quietly
cat("---- Now running all species ----\n")
results <- sapply(species_list, download_phylopic_png, outdir = output_folder, verbose=FALSE)
summary_table <- table(results)
cat("Download summary:\n")
print(summary_table)
