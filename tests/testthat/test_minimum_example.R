library(streambugs)
context("streambugs")


source_test_helpers()


# Note: cannot put model simulation into helper_run_X_example.R because, then it also
#       runs with devtools::document()
if (!exists("minimum.example.setup")) {
    example <- new.env()
    source("run_minimum_example.R", local=example)
    minimum.example.setup <- TRUE
}


test_that("check minimum example feeding links count", {
    n.feeding.links.habitat.expected <- 10
    n.feeding.links.expected <- 3*2*n.feeding.links.habitat.expected
    n.feeding.links <- count.feeding.links(example$model$y.names, example$model$par)
    expect_equal(n.feeding.links.expected, n.feeding.links)
    
    n.habitats <- 5
    n.reaches <- 1
    n.feeding.links.expected <- 1*5*n.feeding.links.habitat.expected
    model <- streambugs.example.model.toy(n.Reaches = n.reaches, n.Habitats = n.habitats)
    n.feeding.links <- count.feeding.links(model$y.names, model$par)
    expect_equal(n.feeding.links.expected, n.feeding.links)
    n.feeding.links <- count.feeding.links(model$y.names$y.names, model$par)
    expect_equal(n.feeding.links.expected, n.feeding.links)
})


test_that("check minimum example numerical solution", {
    .expect_equal_res_matrices("minimum_example_res_R.csv", example$res.R$res)
    .expect_equal_res_matrices("minimum_example_res_C.csv", example$res.C$res)
})


test_that("check minimum example plots", {
    expected.pdf <- "minimum_example_res.pdf"
    .expect_equal_pdf_size_res(example$res.R, expected.pdf, plot.res.func=plot.streambugs)
    .expect_equal_pdf_size_res(example$res.R, expected.pdf)
    .expect_equal_pdf_size_res(example$res.C, expected.pdf)
    expected.pdf <- "minimum_example_foodweb-texts.pdf"
    .expect_equal_pdf_size_foodweb(example$model, expected.pdf, texts=TRUE)
    expected.pdf <- "minimum_example_foodweb-points.pdf"
    .expect_equal_pdf_size_foodweb(example$model, expected.pdf, texts=FALSE)
})


test_that("check minimum example system defintion", {
    expected.dat <- "minimum_example_sysdef.dat"
    .expect_equal_dat_file(example$res.R, expected.dat)
    .expect_equal_dat_file(example$res.C, expected.dat)
})
