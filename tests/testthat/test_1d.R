context("1d clustering")
data(GvHD)
fr <- GvHD[[1]]
trans <- estimateLogicle(fr, c("FL1-H", "FL2-H", "FL3-H", "FL4-H", "FL2-A"))
fr <- flowCore::transform(fr, trans)



test_that("flowClust:SSH-H, 1 mode", {
  chnl <- "SSC-H"
  options("mc.cores" = 1)  
  expect_message(res <- flowClust(fr, varNames = chnl, tol = 1e-7, K = 2), regexp = "serial")
  expect_equal(res@mu, matrix(c(5.27, 6.84), nrow = 2), tol = 5e-2)
  expect_equal(res@w, c(0.286, 0.713), tol = 3e-2)
  
  options("mc.cores" = 4)  
  #relax tol to see the worse fit
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 2)
  expect_equal(res@mu, matrix(c(5.27, 6.84), nrow = 2), tol = 0.3)
  expect_equal(res@w, c(0.286, 0.713), tol = 0.16)
  
  #fiddle with randomStart to change the fit
  res <- flowClust(fr, varNames = chnl, tol = 1e-10, K = 2, randomStart = 10)
  expect_equal(res@mu, matrix(c(5.27, 6.84), nrow = 2), tol = 5e-3)
  expect_equal(res@w, c(0.286, 0.713), tol = 2e-3)
  
  #tmixture
  g <- tmixFilter(parameters = chnl, tol = 1e-10, K = 2, randomStart = 10)
  res.g <- filter(fr, g)
  res.g <- as(res.g, "flowClust")
  res.g@varNames <- as.vector(res.g@varNames)
  expect_equal(res, res.g)
  
  #run through K=1:4
  res <- flowClust(fr, varNames = chnl, tol = 1e-10, K = 1:4, randomStart = 0)
  expect_equal(res[[2]]@mu, matrix(c(5.27, 6.84), nrow = 2), tol = 5e-3)
  expect_equal(res[[2]]@w, c(0.286, 0.713), tol = 6e-4)
  
  
  expect_equal(res[[1]]@mu, matrix(c(6.4), nrow = 1), tol = 6e-3)
  
  #k = 1 is expected to be the best fit
  scores <- sapply(res, slot, "ICL")
  #ICLs decrease as K increases
  expect_true(all(diff(scores) <0))
  
  # par(mfrow=c(1,4))
  # for(obj in res)
  #  hist(obj, fr, main = paste("ICL:", round(obj@ICL),"BIC:", round(obj@BIC)))
  
  
})


test_that("flowClust:FL1-H, 2 mode", {
  
  chnl <- "FL1-H"
  
  res <- flowClust(fr, varNames = chnl, tol = 1e-10, K = 1:4, randomStart = 0)
  
  #K=1
  expect_equal(res[[1]]@mu, matrix(c(0.7057928), nrow = 1), tol = 0.001)
  
  #K=2
  expect_equal(res[[2]]@mu, matrix(c(0.3590909, 1.4315404), nrow = 2), tol = 0.001)
  expect_equal(res[[2]]@w, c(0.6837035, 0.3162965), tol = 7e-4)
  
  #K=3
  expect_equal(res[[3]]@mu, matrix(c(0.001814324, 0.515071157, 1.424679742), nrow = 3), tol = 0.001)
  expect_equal(res[[3]]@w, c(0.2218212, 0.4538276, 0.3243512), tol = 7e-4)
  
  
  #2 mode is best fit
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_gt(scores.diff[1], 0) #1st is pos
  expect_true(all(scores.diff[-1] <0))# the rest are neg
  
  par(mfrow=c(1,4))
  for(obj in res)
   hist(obj, fr, main = paste("ICL:", round(obj@ICL),"BIC:", round(obj@BIC)))
  
})

test_that("flowClust:FL2-H, 3 mode", {
  chnl <- "FL2-H"
  
  res <- flowClust(fr, varNames = chnl, tol = 1e-10, K = 1:4, randomStart = 0)
  expect_equal(res[[1]]@mu, matrix(c(1.820132), nrow = 1), tol = 0.001)
  
  expect_equal(res[[2]]@mu, matrix(c(1.610172, 3.148424), nrow = 2), tol = 0.001)
  expect_equal(res[[2]]@w, c(0.769645, 0.230355), tol = 7e-4)
  
  expect_equal(res[[3]]@mu, matrix(c(-0.09441319, 1.21476765, 1.91313389), nrow = 3), tol = 0.001)
  expect_equal(res[[3]]@w, c(0.1041290, 0.6273802, 0.2684908), tol = 7e-4)
  
  
  
  #3 mode is best fit
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_gt(scores.diff[2], 0) #2nd is pos
  expect_true(all(scores.diff[-2] < 0))# the rest are neg
  # par(mfrow=c(1,4))
  # for(obj in res)
  #  hist(obj, fr, main = paste("ICL:", round(obj@ICL),"BIC:", round(obj@BIC)))
})

test_that("flowClust:FL2-A, 3 mode", {
  # chnl <- "FL2-A"
  
  # res <- flowClust(fr, varNames = chnl, tol = 1e-10, K = 1:4, randomStart = 0)
  
  #3 mode is best fit
  # scores <- sapply(res, slot, "ICL")
  # scores.diff <- diff(scores)
  # expect_gt(scores.diff[2], 0) #2nd is pos
  # expect_true(all(scores.diff[-2] < 0))# the rest are neg
  # par(mfrow=c(1,4))
  # for(obj in res)
  #  hist(obj, fr, main = paste("ICL:", round(obj@ICL),"BIC:", round(obj@BIC)))
})


