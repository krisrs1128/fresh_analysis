
get_data_dir <- function(){
  wd <- getwd()
  if (str_detect(wd, "laurasymul"))
    data_dir <- "/Users/laurasymul/Dropbox/Academia/Projects/2024 FRESH Multiomics R01/data_sets/"
  else if (str_detect(wd, "krissankaran"))
    data_dir <- "/Users/krissankaran/Desktop/collaborations/fresh_surveys/data/"
  else if (str_detect(wd, "vastm"))
    data_dir <- "/Users/vastm/OneDrive - UCL/Documents/Thesis/Shared Madeline - Laura/02 Projects/FRESH HIV/data_sets/"
  else if (str_detect(wd, "camst"))
    data_dir <- "C:/Users/camst/Dropbox/data_sets/"
  else if (str_detect(wd, "vanbenedena"))
    data_dir <- "/Users/vanbenedena/Dropbox/2024 FRESH Multiomics R01/data_sets/"
  else if (str_detect(wd, "arnalot"))
    data_dir <- "/Users/arnalot/Dropbox/2024 FRESH Multiomics R01/data_sets/"
  else
    stop("get_data_dir() not specified for this user/computer")
  data_dir
}

get_redcap_dir <- function(){
  str_c(get_data_dir(), "REDCap exports/raw_export_20240730/")
}
