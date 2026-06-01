# Build and save the toy example MAE.
#   Rscript src/preprocess/build_example_mae.R
library(here)
library(fs)
r <- path(here())

source(r / "src/data/demo/simulate.R")
mae <- build_example_mae(seed = 202605)
dir.create("data", showWarnings = FALSE)
saveRDS(mae, r / "data/example_mae.rds")

cat("Saved data/example_mae.rds\n")
print(mae)