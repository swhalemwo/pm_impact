library(data.table)

# Load data and source logic
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Run matching in full return mode
# dt has multiple years, let's just use the whole thing
print("Starting matching (full return mode for whole dt)...")
start_time <- Sys.time()
full_matches <- match_treatment_history(dt, L = 3, history_col = "exn_pm_cum", return_summary = FALSE)
end_time <- Sys.time()

print(paste("Full matching completed in:", end_time - start_time))
print(paste("Number of matches found:", nrow(full_matches)))
print(head(full_matches))
