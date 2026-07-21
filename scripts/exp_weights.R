

dt2 <- data.table(stratum = c(1,1,1,1,2,2,2,2,2,2),
                  treat1 = c(0,0,1,1,0,1,1,0,0,1),
                  treat2 = c(1,0,1,0,0,0,1,0,1,0))

dt_s2 <- dt2[, .(msT1 = sum(treat1), msT2 = sum(treat2), msC = sum(treat1 == 0 & treat2 == 0)), stratum]

## need global ratios: have 3 groups -> 3 ratios (2+1), for N treatments it's N + N-1 ratios
rat_c_t1 = sum(dt_s2$msC)/sum(dt_s2$msT1)
rat_c_t2 = sum(dt_s2$msC)/sum(dt_s2$msT2)
rat_t1_t2 = sum(dt_s2$msT1)/sum(dt_s2$msT2)


treat_cols <- c("treat1", "treat2")
n_tk   <- sapply(treat_cols, \(tc) sum(dt2[[tc]]))          ## per-treatment counts
n_notk <- sapply(treat_cols, \(tc) sum(dt2[[tc]] == 0))     ## complements
R      <- n_notk / n_tk                                      ## target: (not treat_k)/(treat_k)
p      <- R / (1 + R)
print(p)




cem_multitreat_weights <- function(dt_s, p_vec) {
  K <- length(p_vec)
  n <- nrow(dt_s)
  A <- matrix(0, K, n)
  for (k in seq_len(K))
    A[k, ] <- dt_s[[names(p_vec)[k]]] - p_vec[k]
  AAt_inv <- solve(A %*% t(A))
  w <- rep(1, n) - crossprod(A, AAt_inv %*% (A %*% rep(1, n)))
  as.numeric(w)
}


dt2[, w := cem_multitreat_weights(.SD, p), stratum]

## my check
dt2[, .(w1 = sum(w * (treat1 ==1))/sum(treat1==1),
        w2 = sum(w * (treat2 ==1))/sum(treat2==1),
        wc = sum(w *(treat1 == 0 & treat2 == 0))/sum(treat1 == 0 & treat2 == 0)), stratum] %>%
.[, `:=`(r1 = w1/wc, r2= w2/wc)]


## AI check
dt2[, .(r1 = sum(w * treat1) / sum(w * (1 - treat1)),
        r2 = sum(w * treat2) / sum(w * (1 - treat2))), by = stratum]
