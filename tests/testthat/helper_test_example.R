

#' Read reference CSV simulation results and compare with the simulated ones.
.expect_equal_res_matrices <- function(csv.expected, res.matrix) {
  res.expected = as.matrix(read.csv(csv.expected, check.names=FALSE))
  testthat::expect_equal(dim(res.expected),dim(res.matrix))
  testthat::expect_equal(colnames(res.expected), colnames(res.matrix))
  # compare values with rel. tolerance of 1e-6
  vec.expected = as.vector(unlist(res.expected))
  vec.matrix = as.vector(unlist(res.matrix))
  vec.rel.err = mapply(function(x,y) if (x == 0) abs(y) else abs((x-y)/x), vec.expected, vec.matrix)
  vec.zero = rep(0,prod(dim(res.expected)))
  testthat::expect_equal(vec.rel.err, vec.zero, tolerance=1e-6)
}


.streambugs.sys.def.warning.handler <- function(w) {
  if(any(grepl("no (envcond|parameters) for .*", w))) {
    invokeRestart("muffleWarning")
  }
}


#' Plot streambugs results to PDF and test if PDF has same size as the reference.
.expect_equal_pdf_size_res <- function(res.list, expected.pdf, plot.res.func=plot) {
  res.pdf <- "temp_example_res.pdf"
  pdf(res.pdf, width=8, height=6)
  withCallingHandlers(
    plot.res.func(res.list$res, res.list$args$par, res.list$args$inp),
    warning = .streambugs.sys.def.warning.handler
  )
  dev.off()
  .expect_equal_pdf_size(res.pdf, expected.pdf)
}


#' Plot model foodweb to PDF and test if PDF has same size as the reference.
.expect_equal_pdf_size_foodweb <- function(model.list, expected.pdf, ...) {
  temp.pdf <- "temp_example_foodweb.pdf"
  streambugs::foodweb.plot(model.list$y.names, model.list$par, file=temp.pdf, ...)
  .expect_equal_pdf_size(temp.pdf, expected.pdf)
}


#' Test if temporary PDF has same size as the reference PDF with a small relative
#' tolerance for an OS-dependent PDF size variation. Cleanup temporary PDF afterwards.
.expect_equal_pdf_size <- function(temp.pdf, ref.pdf, reltol=0.005) {
  tryCatch({
    ref.pdf.size <- file.info(ref.pdf)$size
    temp.pdf.size <- file.info(temp.pdf)$size
    abstol = reltol*ref.pdf.size
    testthat::expect_equal(ref.pdf.size, temp.pdf.size, tolerance=abstol)
  },
  finally={
    file.remove(temp.pdf)
  })
}


.read.file.to.char <- function(file.name) {
  # standardize EOLs: Windows to Unix
  gsub('\r\n', '\n', readChar(file.name, file.info(file.name)$size))
}


.expect_equal_dat_file <- function(res.list, expected.sysdef.dat) {
  # streambugs.get.sys.def
  sys.def <- withCallingHandlers(
    streambugs::streambugs.get.sys.def(
        y.names=res.list$args$y.names,
        par=res.list$args$par,
        inp=res.list$args$inp),
    warning = .streambugs.sys.def.warning.handler
  )
  # streambugs.write.sys.def
  res.dat <- "temp_example_sysdef.dat"
  streambugs::streambugs.write.sys.def(sys.def, file=res.dat)
  # read both files to string and compare
  expected.dat.char <- .read.file.to.char(expected.sysdef.dat)
  tryCatch({
    testthat::expect_equal(expected.dat.char, .read.file.to.char(res.dat))
  },
  finally={
    file.remove(res.dat)
  })
}
