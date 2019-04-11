context("2d clustering")
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
fs <- read.flowSet(fcsFiles)
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "warn")


test_that("flowClust 2d: prior", {
  set.seed(1)
  require(openCyto)
  
  prior <- prior_flowClust(fs ,c("FSC-A", "SSC-A"), K = 2)
  options("mc.cores" = 1)  
  g <-
  	flowClust.2d(
  		fs[[1]],
  		"FSC-A",
  		"SSC-A",
  		K = 2,
  		target = c(1e5, 5e4),
  		prior = prior,
  		usePrior = "yes"
  	)
  
  expect_equal(g@cov, structure(c(1173858741.42853, 523766402.261612, 523766402.261612, 625519945.097251)
                                , .Dim = c(2L, 2L)
                                , .Dimnames = list(c("FSC-A", "SSC-A"), c("FSC-A", "SSC-A")))
              )
  expect_equal(g@distance, 2.94, tolerance =1e-3)
  expect_equal(g@mean, c("FSC-A" = 87638.359375, "SSC-A" = 62425.46875))
})

