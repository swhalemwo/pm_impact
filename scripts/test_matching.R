library(data.table)
library(pryr)
source("scripts/mymatch.R")

# Load small subset for testing
# dt_test is already in the global env (or we assume it's available)
# dt_test <- ... (ensure it is loaded)

run_test_match <- function(n_treated) {
  message("--- Starting test with n_treated = ", n_treated, " ---")
  gc()
  mem_start <- mem_used()
  
  # Sample a subset of treated units
  treated_ids <- unique(dt_test[exn_pm > 0, PeopleID])
  sample_ids <- sample(treated_ids, min(n_treated, length(treated_ids)))
  dt_sample <- dt_test[PeopleID %in% sample_ids]
  
  # Profile match
  mem_pre_match <- mem_used()
  match_res <- match_treatment_history(dt_sample, L = 3)
  mem_post_match <- mem_used()
  
  mem_delta <- mem_post_match - mem_pre_match
  message("Memory used by match: ", mem_delta)
  
  return(list(delta = mem_delta, results = match_res))
}

# Iterative execution
# Start with n=10, then double
sizes <- c(10, 20, 40, 80, 160, 320, 640, 1280, 2560)
for (n in sizes) {
  res <- run_test_match(n)
  if (res$delta > 500 * 1024^2) { # Halt if usage exceeds 500MB delta
     stop("Memory usage exceeded threshold!")
  }
}
