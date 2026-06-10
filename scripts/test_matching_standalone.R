library(data.table)

# Load data and source logic
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Run matching in summary mode
print("Starting matching (summary mode)...")
start_time <- Sys.time()
match_summary <- match_treatment_history(dt, L = 3, history_col = "exn_pm_cum", return_summary = TRUE)
end_time <- Sys.time()

print(match_summary)
print(paste("Summary completed in:", end_time - start_time))
