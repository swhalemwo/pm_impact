library(data.table)

# Load data and source logic
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# 1. Run Summary Mode
print("Running Summary Mode...")
summary_res <- gd_match_history(dt, L = 3, return_summary = TRUE)
print(summary_res)

# 2. Run Full Match Mode
print("Running Full Match Mode...")
full_matches <- gd_match_history(dt, L = 3, return_summary = FALSE)

# 3. Calculate Summary from Full Matches
full_summary <- data.table(
  total_treated = uniqueN(dt[exn_pm > 0, PeopleID]), # Adjusted for exn_pm check
  treated_matched = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = mean(full_matches[, .N, by = treated_id]$N),
  total_matches = nrow(full_matches)
)

print("Summary from Full Matches:")
print(full_summary)

# Compare
if(all.equal(summary_res$total_matches, full_summary$total_matches) &&
   all.equal(summary_res$treated_matched, full_summary$treated_matched)) {
  print("SUCCESS: Summary and Full Match counts match.")
} else {
  stop("FAILURE: Counts do not match.")
}
