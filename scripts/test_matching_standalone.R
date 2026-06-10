library(data.table)

# Load data and source logic
dt <- fread("~/Dropbox/phd/papers/foster/dt_test.csv")
source("~/Dropbox/phd/papers/foster/scripts/mymatch.R")

# Run matching
print("Starting matching...")
start_time <- Sys.time()
matched <- match_treatment_history(dt, L = 3)
end_time <- Sys.time()

# Validation
print(paste("Matching completed in:", end_time - start_time))
print(paste("Number of matches found:", nrow(matched)))

# Basic correctness check
# Correct indexing for matched data.table
set.seed(42)
sample_rows <- matched[sample(.N, 10)]

check_match <- function(m_row) {
  t_id <- m_row$treated_id
  c_id <- m_row$control_id
  yr <- m_row$year
  
  # Get histories for both units in the year leading up to treatment
  # History for L=3: years y-3, y-2, y-1
  get_hist <- function(id, y) {
    # Sort ensure order is correct for paste
    dt[PeopleID == id & year %in% (y-3):(y-1)][order(year), as.integer(exn_pm > 0)]
  }
  
  h1 <- get_hist(t_id, yr)
  h2 <- get_hist(c_id, yr)
  
  # Check if histories are valid length and equal
  if(length(h1) < 3 || length(h2) < 3) return(FALSE)
  return(all(h1 == h2))
}

checks <- sapply(1:nrow(sample_rows), function(i) check_match(sample_rows[i]))

if (all(checks == TRUE)) {
  print("Validation passed: All sampled matches have identical treatment history.")
} else {
  stop("Validation failed: Mismatched treatment history found.")
}
