context("1d clustering on gated data")
library(flowWorkspace)
dataDir <- system.file("extdata",package="flowWorkspaceData")
gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
Sys.setenv("_R_CHECK_LIMIT_CORES_"="warn")
options("mc.cores" = 4)
test_that("singlets node:", {
  fr <- getData(gs[[1]], "singlets")
  pd <- pData(parameters(fr))
  pd <- pd[!is.na(pd[["desc"]]), ]
  channels <- pd[["name"]]
  fr <- fr[, channels]
  
  chnl <- "<B710-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0,)
  
  #3 mode and expect k = 3 is best
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(all( scores.diff[1:2] > 0))
  expect_true(all( scores.diff[3] <0))
  
  chnl <- "<R660-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #1 mode, no good separation
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(all( scores.diff <0))
  
  chnl <- "<R780-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #2 mode
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(scores.diff[1] > 0)
  expect_true(which.max(scores.diff) == 1)
  
  chnl <- "<V450-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #2 mode
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(scores.diff[1] > 0)
  expect_true(which.max(scores.diff) == 1)
  cd3.score <- scores.diff[1]
  
  chnl <- "<V545-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #2 mode
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(scores.diff[1] > 0)
  expect_true(which.max(scores.diff) == 1)
  HLA.score <- scores.diff[1]
  
  chnl <- "<G560-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #1 mode
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(all(scores.diff < 0))
  
  chnl <- "<G780-A>"
  res <- flowClust(fr, varNames = chnl, tol = 1e-5, K = 1:4, randomStart = 0, min.count = -1, max.count = -1, trans = 0)
  #2 mode
  scores <- sapply(res, slot, "ICL")
  scores.diff <- diff(scores)
  expect_true(scores.diff[1] > 0)
  expect_true(which.max(scores.diff) == 1)
  cd45.score <- scores.diff[1]
  
  # par(mfrow=c(1,4))
  # for(obj in res)
  #  hist(obj, fr, main = paste("ICL:", round(obj@ICL),"BIC:", round(obj@BIC)))
  # 
  # 
})

