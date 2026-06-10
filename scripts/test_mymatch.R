library(data.table)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Load data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")

# 1. Test Summary Mode
print("--- Testing Summary Mode ---")
summary_res <- match_treatment_history(dt, L = 3, return_type = "summary")
print(summary_res)

# 2. Test Matches Mode
print("--- Testing Matches Mode ---")
matches <- match_treatment_history(dt, L = 3, return_type = "matches")
if(all(c("control_id", "treated_id", "year") %in% names(matches)) && 
   !any(c("match_key", "match_tier") %in% names(matches))) {
  print("SUCCESS: Matches mode returns only base columns.")
} else {
  stop("FAILURE: Matches mode returns incorrect columns.")
}

# 3. Test Matches-Info Mode
print("--- Testing Matches-Info Mode ---")
matches_info <- match_treatment_history(dt, L = 3, return_type = "matches-info")
if(all(c("match_key", "match_tier") %in% names(matches_info))) {
  print("SUCCESS: Matches-Info returns requested debugging columns.")
} else {
  stop("FAILURE: Missing debugging columns.")
}

# 4. Check Tier counts
print("Tier distribution:")
print(table(matches_info$match_tier))

print("All tests passed.")
