library(streambugs)
context("streambugs")


source_test_helpers()


# Note: cannot put model simulation into helper_run_X_example.R because, then it also
#       runs with devtools::document()
if (!exists("extended.example.setup")) {
    example <- new.env()
    source("run_extended_example.R", local=example)
    extended.example.setup <- TRUE
}


test_that("check extended example feeding links count", {
    n.feeding.links.habitat.expected <- 290
    n.feeding.links.expected <- 1*6*n.feeding.links.habitat.expected
    n.feeding.links <- count.feeding.links(example$model$y.names, example$model$par)
    expect_equal(n.feeding.links.expected, n.feeding.links)
})


test_that("check extended example numerical solution", {
    .expect_equal_res_matrices("extended_example_res_C.csv", example$res.C$res)
})


test_that("check extended example plots", {
    expected.pdf <- "extended_example_res.pdf"
    .expect_equal_pdf_size_res(example$res.C, expected.pdf)
    expected.pdf <- "extended_example_foodweb.pdf"
    .expect_equal_pdf_size_foodweb(example$model, expected.pdf)
})


test_that("check extended example system defintion", {
    expected.dat <- "extended_example_sysdef.dat"
    .expect_equal_dat_file(example$res.C, expected.dat)
})
