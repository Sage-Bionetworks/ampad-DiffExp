f.pvalue <- function (dat, mod, mod0, tol = 1e-18) 
{
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0, m)
  Id <- diag(n)
  resid <- dat %*% (Id - mod %*% solve(t(mod) %*% mod, tol = tol) %*% 
                      t(mod))
  rss1 <- rowSums(resid * resid)
  rm(resid)
  resid0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0, tol = tol) %*% 
                       t(mod0))
  rss0 <- rowSums(resid0 * resid0)
  rm(resid0)
  fstats <- ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
  p <- 1 - pf(fstats, df1 = (df1 - df0), df2 = (n - df1))
  return(p)
}