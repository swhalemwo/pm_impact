library(data.table)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Load data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")

# 1. Test Summary Mode
print("--- Testing Summary Mode ---")
summary_res <- match_treatment_history(dt, L = 3, return_summary = TRUE)
print(summary_res)

# 2. Test Full Match Mode
print("--- Testing Full Match Mode ---")
full_matches <- match_treatment_history(dt, L = 3, return_summary = FALSE)

# 3. Comprehensive Consistency Check
print("--- Running Rigorous Consistency Check ---")
# avg_matches_per_treated for summary mode is calculated as:
# total_matches / uniqueN(treated_set)
# For full_matches, it is: mean(matches$N) from the internal calculation
# These are identical.

full_summary <- data.table(
  total_treated = uniqueN(dt[exn_pm > 0, PeopleID]),
  treated_matched = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = nrow(full_matches) / uniqueN(full_matches$treated_id),
  total_matches = nrow(full_matches)
)

print(full_summary)

# Compare metrics using identical logic
stopifnot(all.equal(as.numeric(summary_res$total_matches), as.numeric(full_summary$total_matches)))
stopifnot(all.equal(as.numeric(summary_res$treated_matched), as.numeric(full_summary$treated_matched)))
stopifnot(all.equal(as.numeric(summary_res$avg_matches_per_treated), as.numeric(full_summary$avg_matches_per_treated)))

print("SUCCESS: Summary and Full Match metrics match exactly.")
print("All tests passed.")
