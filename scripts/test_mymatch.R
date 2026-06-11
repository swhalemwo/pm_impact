library(data.table)
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Load data
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")

# 1. Test Matching (Summary Mode)
print("--- Testing Summary Mode ---")
summary_res <- match_treatment_history(dt, L = 3, return_type = "summary")
print(summary_res)

# 2. Test Full Match Mode
print("--- Testing Full Match Mode ---")
full_matches <- match_treatment_history(dt, L = 3, return_type = "matches")
print(paste("Full matches rows:", nrow(full_matches)))


if(all(c("control_id", "treated_id", "year") %in% names(full_matches)) && 
   !any(c("match_key", "match_tier") %in% names(full_matches))) {
  print("SUCCESS: Matches mode returns only base columns.")
} else {
  stop("FAILURE: Matches mode returns incorrect columns.")
}

# 3. Consistency Check (Replicate Summary from Full)
print("--- Running Rigorous Consistency Check ---")
full_summary <- data.table(
  total_treated = uniqueN(dt[exn_pm > 0, PeopleID]),
  treated_matched = uniqueN(full_matches$treated_id),
  avg_matches_per_treated = mean(full_matches[, .N, by = treated_id]$N),
  total_matches = nrow(full_matches)
)
stopifnot(all.equal(as.numeric(summary_res$total_matches), as.numeric(full_summary$total_matches)))
stopifnot(all.equal(as.numeric(summary_res$treated_matched), as.numeric(full_summary$treated_matched)))
print("SUCCESS: Consistency check passed.")

# 4. Test Matches-Info Mode
print("--- Testing Matches-Info Mode ---")
matches_info <- match_treatment_history(dt, L = 3, return_type = "matches-info")
stopifnot(all(c("match_key", "match_tier") %in% names(matches_info)))
print("SUCCESS: Matches-Info returns debugging columns.")

# 5. Test Refinement (CEM) on Full Data
print("--- Testing Refinement (CEM) on Full Data ---")
# Use the full match output
dt_refined <- refine_matches(full_matches, dt, 
                             covariates = c("BirthYear", "exn_all_cum"),
                             cutpoints = list(BirthYear = 8, exn_all_cum = 8))

if("subclass" %in% names(dt_refined)) {
  print("SUCCESS: Refinement returned subclass assignments.")
  print(head(dt_refined))
} else {
  stop("FAILURE: Refinement failed.")
}

# 4. Test number of tiers
dt_tiercnts <- matches_info[, .(nbr_tiers = uniqueN(match_tier)), treated_id]
if (max(dt_tiercnts$nbr_tiers) > 1) {
    stop("some treated are in multiple tiers")
} else {
    print("SUCCESS: all treated are in only one tier")
}
print("All tests passed.")
