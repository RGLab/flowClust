context("2d clustering")
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
fs <- read.flowSet(fcsFiles)
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "warn")


test_that("flowClust 2d: prior", {
  RNGversion("3.5")#quick fix to pass the test against the results that were previously generated from older R (<=2018/11/06)
  set.seed(1)
  require(openCyto)
  
  prior <- prior_flowclust(fs ,c("FSC-A", "SSC-A"), K = 2)
  options("mc.cores" = 1)  
  g <-
  	gate_flowclust_2d(
  		fs[[1]],
  		"FSC-A",
  		"SSC-A",
  		K = 2,
  		target = c(1e5, 5e4),
  		prior = prior,
  		usePrior = "yes"
  	)
  
  expect_equal(g@cov, structure(c(1173858741, 523766402, 523766402, 625519945)
                                , .Dim = c(2L, 2L)
                                , .Dimnames = list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A")))
               , tolerance = 1e-6
               )
  expect_equal(g@distance, 2.94, tolerance =1e-3)
  expect_equal(g@mean, c("FSC-A" = 87638, "SSC-A" = 62425), tolerance = 1e-5)
})

