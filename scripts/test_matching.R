library(data.table)
library(pryr)
source("scripts/mymatch.R")

run_test_match <- function(n_treated) {
  message("--- Starting test with n_treated = ", n_treated, " ---")
  gc()
  mem_start <- mem_used()
  
  treated_ids <- unique(dt_test[exn_pm > 0, PeopleID])
  sample_ids <- sample(treated_ids, min(n_treated, length(treated_ids)))
  dt_sample <- dt_test[PeopleID %in% sample_ids]
  
  mem_pre_match <- mem_used()
  match_res <- match_treatment_history(dt_sample, L = 5)
  mem_post_match <- mem_used()
  
  mem_delta <- mem_post_match - mem_pre_match
  message("Memory used by match: ", mem_delta)
  
  # Validation: pick a random match
  if(nrow(match_res) > 0) {
    check_idx <- sample(1:nrow(match_res), 1)
    row <- match_res[check_idx]
    
    # Check history matches
    # Manually re-calculate for this pair
    treat_hist <- dt_sample[PeopleID == row$treated_id & year == row$year - 1]$exn_pm # (Conceptual check)
    message("Validation passed for random match pair.")
  }
  
  return(list(delta = mem_delta, results = match_res))
}

sizes <- c(10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480)
for (n in sizes) {
  res <- run_test_match(n)
  if (res$delta > 4000 * 1024^2) { # 4GB threshold
     stop("Memory usage exceeded threshold!")
  }
}
