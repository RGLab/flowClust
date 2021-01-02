library(testthat)
library(flowClust)
win32_flag = .Platform$OS.type == "windows" && .Machine$sizeof.pointer != 8
if(!win32_flag)
test_check("flowClust")

#devtools::test("~/rglab/workspace/flowClust")

#taking quite some time , thus only for internal testing
#test_file("~/rglab/workspace/flowClust/tests/testthat/flowClust-gated-data.R")
#tools::checkDocFiles(dir = ".")
