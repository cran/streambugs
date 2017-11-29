# README

Developers instructions for building, installing and testing Streambugs R
package sources.

## Installation

### Prerequisites

You need a C compiler (recommended GCC >= 5 or LLVM >= 7); cf.

    $ gcc --version

You also need an R (>= 3.0.2) installation on your computer; run:

    $ R
    > install.packages(c('testthat', 'devtools', 'roxygen2', 'deSolve'), dependencies=TRUE)

which should download, build and install all R tools needed to build C/C++ based
packages.

> #### Attention
> installing above R packages, together with all suggested dependencies (the
> `dependencies=TRUE` keword argument), takes some time, especially for the
> `devtools` package. You might opt for not installing suggested dependencies,
> and install the missing ones as you discover them.

### Build, install and test

Run `R` shell and load `devtools` library:

    $ R
    > library(devtools)

To create docs, build package and install it run in this folder:

    > devtools::document()
    > pkg_src = devtools::build()
    > install.packages(pkg_src, repos=NULL, type="source")

#### Run tests

To run tests, with the package as a working directory, run:

    > devtools::test()

#### Check for CRAN submission

To check for CRAN submission, with the package as a working directory, run:

    > devtools::check()

##### Check archive for submission

Build will create a zipped TAR file `../streambugs_PKGVER.tar.gz`,
where `PKGVER` is the current package version. You can check this package for
CRAN submission by running in shell:

    $ R CMD check --as-cran ../streambugs_$(VERSION).tar.gz

This creates the `streambugs.Rcheck/` folder which contains also a PDF file
`streambugs-manual.pdf`, that is the R vignette of the package.

Equivalently, from `R` shell, after the build, run:

    > devtools::check(pkg_src)

### Cleanup

Remove build, packaging files, and CRAN check folder, i.e., relative to this
folder, following files:

    src/*.o
    src/*.so
    man/*.Rd
    ../streambugs_0.1.tar.gz
    streambugs.Rcheck/

To uninstall the package from your system run in `R` shell:

    > remove.packages("streambugs")
