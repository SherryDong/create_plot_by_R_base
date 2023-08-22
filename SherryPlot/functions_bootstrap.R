# general bootstrap for any statistics and comparison (paired t-test)
bootstrap_stat <- function(
  df, n_boot, func, score_names,
  compare_T = F, compare_cols = NULL,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  need_se = F,
  seed = 0
){
  set.seed(seed)
  res_boot <- matrix(nrow = n_boot, ncol = length(score_names))
  bootstrap_list <- list()
  
  for (k in 1:n_boot){
    idx_k <- sample.int(nrow(df), replace = T)
    df_k <- df[idx_k, ]
    res_k <- func(df_k)
    res_boot[k, ] <- res_k
  }
  bootstrap_list$bootstrap_value <- res_boot
  
  stat_boot <- apply(res_boot, 2, function(x){quantile(x, probs = probs, na.rm = T)})
  rownames(stat_boot) <- lapply(probs, function(x){sprintf('ci_%.3f', x)})
  bootstrap_list$stat_boot <- stat_boot
  
  if (compare_T){
    xy <- res_boot[, compare_cols]
    res_ttest <- t.test(xy[, 1], xy[, 2], paired = T)
    bootstrap_list$compare_T <- c(
      p.value=res_ttest$p.value[[1]],
      estimate=res_ttest$estimate[[1]]
    )
    
  }
  if (need_se){
    se_boot <- apply(res_boot, 2, sd)
    bootstrap_list$stat_boot <- cbind(stat_boot, se=se_boot)
  }
  return(bootstrap_list)
}