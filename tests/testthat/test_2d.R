context("2d clustering")
flowDataPath <- system.file("extdata", package = "flowWorkspaceData")
fcsFiles <- list.files(pattern = "CytoTrol", flowDataPath, full = TRUE)
fs <- read.flowSet(fcsFiles)
Sys.setenv("_R_CHECK_LIMIT_CORES_" = "warn")


test_that("flowClust 2d: prior", {
  
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
  
  # expect_equal(g@boundaries, matrix(c(5.27, 6.84), nrow = 2))
})

