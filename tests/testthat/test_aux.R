context("streambugs")


test_that("sparsify.inpmatch", {
    m <- array(
        data = c(1, NA, 3, NA, 5, 6),
        dim = c(3, 2),
        dimnames = list(c("a", "b", "c"), NULL))
    expected.sm <- rbind(
        c(1, 1, 1),
        c(3, 3, 1),
        c(5, 2, 2),
        c(6, 3, 2))
    colnames(expected.sm) <- c("inpind", "i", "j")
    rownames(expected.sm) <- rownames(m)[expected.sm[,"i"]]
    actual.sm <- streambugs:::sparsify.inpmatch(m)
    testthat::expect_equal(actual.sm, expected.sm)
})


test_that("check.inpind.parval.required", {
    test.data = list(
        d2 = list(
            parvals = matrix(
                data = c(1:2, NA, 4:6), nrow = 2, ncol = 3,
                dimnames = list(c("x1", "x2"), c("p1", "p2", "p3"))),
            inpinds = matrix(
                data = c(3, 1, 2), nrow = 1, ncol = 3,
                dimnames = list(NULL, c("inpinds", "i", "j")))
        ),
        d1 = list(
            parvals = stats::setNames(c(1, NA, 3), c("p1", "p2", "p3")),
            inpinds = matrix(
                data = c(2, 2), nrow = 1, ncol = 2,
                dimnames = list(NULL, c("inpinds", "i")))

        )
    )
    for (inpinds.parvals in test.data) {
        # all good, even though NA for var #1 "p2" because it's in input pars instead
        testthat::expect_silent(
            streambugs:::check.inpind.parval.required(
                inpinds.parvals = inpinds.parvals,
                required = c("p2", "p1"))
        )
        # missing "p2" value for var #1
        testthat::expect_warning(
            streambugs:::check.inpind.parval.required(
                inpinds.parvals = list(
                    inpinds = inpinds.parvals$inpinds[NULL,],
                    parvals = inpinds.parvals$parvals),
                required = c("p2"))
        )
        # required "p4" is not in prarameters names
        testthat::expect_warning(
            streambugs:::check.inpind.parval.required(
                inpinds.parvals = inpinds.parvals,
                required = c("p4"))
        )
    }
})

.expect_equal.inpind.parval <- function(
    expected.fn.prefix,
    actual.list,
    model,
    model.name,
    par.names,
    par.names.required=NA,
    par.defaults=NA)
{
    # logic tests
    # TODO: test actual.list$inpind names
    # colnames(actual.list$parvals)[unique(actual.list$inpinds[,"j"])] vs. model$inp "_par" suffix (use regexp)
    # rownames(actual.list$parvals)[actual.list$inpinds[,"i"]] vs. model$inp w/o "_par" suffix (use regexp)
    testthat::expect_equal(colnames(actual.list$parvals), par.names)
    if (!all(is.na(par.names.required)))
        testthat::expect_true(all(!is.na(actual.list$parvals[,par.names.required])))
    # OPT: test par.defaults values for pars not in model$par

    # regression test
    expected.fn <- paste0(expected.fn.prefix, "-", model.name, ".rds")
    expected.list <- readRDS(expected.fn)
    testthat::expect_equal(actual.list, expected.list)
}


models <- list(
    # Note: at least 6 reaches to have 2 input list entires ($Reach3_w and $Reach6_w)
    minimum_example = streambugs::streambugs.example.model.toy(n.Reaches = 6,
                                                               n.Habitats = 2),
    extended_example = streambugs::streambugs.example.model.extended())


test_that("get.inpind.parval.envcond.reach", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        # add multiple inputs for a single var
        for (inp.name in names(model$inp)) {
            model$inp[[sub("_w", "_L", inp.name)]] <- model$inp[[inp.name]]
        }
        for (par.names in list(c("w","L"), c("w"))) {
            par.names.required = par.names
            par.envcond.reach <- streambugs:::get.inpind.parval.envcond.reach(
                par.names = par.names,
                y.names   = model$y.names,
                par       = model$par,
                inp       = model$inp,
                required  = par.names.required)
            .expect_equal.inpind.parval(
                paste0("get_inpind_parval_envcond_reach-", length(par.names), "par"),
                par.envcond.reach,
                model,
                model.name,
                par.names,
                par.names.required = par.names.required)
        }
    }
})


test_that("get.inpind.parval.envcond.habitat", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        par.names <- c("T",
                       "I0",
                       "fshade",
                       "CP",
                       "CN",
                       "DSusPOM",
                       "tau",
                       "taucrit",
                       "tempmaxK",
                       "currentms",
                       "orgmicropollTU",
                       "saprowqclass",
                       "fA",
                       "DFish")
        par.names.required <- c("T")
        par.defaults <- c(fA=1,DFish=0)
        par.envcond.habitat <- streambugs:::get.inpind.parval.envcond.habitat(
            par.names = par.names,
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            required  = par.names.required,
            defaults  = par.defaults)
        .expect_equal.inpind.parval(
            "get_inpind_parval_envcond_habitat",
            par.envcond.habitat,
            model,
            model.name,
            par.names,
            par.names.required = par.names.required,
            par.defaults = par.defaults)
    }
})


test_that("get.inpind.parval.initcondinput", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        # Note: quering for "w" likely doesn't make sense semantically, but it works and
        #       does actually test the $inpinds part in the minimum example case
        for (par.defaults in list(c(Dini=0, Input=0, w=0), c(Dini=0), c(Input=0))) {
            par.names = names(par.defaults)
            par.initcondinput <- streambugs:::get.inpind.parval.initcondinput(
                par.names = par.names,
                y.names   = model$y.names,
                par       = model$par,
                inp       = model$inp,
                defaults  = par.defaults)
            infix <- paste(par.names, collapse = "_")
            .expect_equal.inpind.parval(
                paste0("get_inpind_parval_initcondinput-", infix),
                par.initcondinput,
                model,
                model.name,
                par.names,
                par.defaults = par.defaults)
        }
    }
})

test_that("get.inpind.parval.initcondinput", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        # Note: quering for "w" likely doesn't make sense semantically, but it works and
        #       does actually test the $inpinds part in the minimum example case
        for (par.defaults in list(c(Dini=0, Input=0, w=0), c(Dini=0), c(Input=0))) {
            par.names = names(par.defaults)
            par.initcondinput <- streambugs:::get.inpind.parval.initcondinput(
                par.names = par.names,
                y.names   = model$y.names,
                par       = model$par,
                inp       = model$inp,
                defaults  = par.defaults)
            infix <- paste(par.names, collapse = "_")
            .expect_equal.inpind.parval(
                paste0("get_inpind_parval_initcondinput-", infix),
                par.initcondinput,
                model,
                model.name,
                par.names,
                par.defaults = par.defaults)
        }
    }
})


test_that("get.single.par.inpind.parval", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        par.initcond.expected <- streambugs:::get.inpind.parval.initcondinput(
            par.names = "Dini",
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            defaults  = c(Dini=0))
        par.input.expected <- streambugs:::get.inpind.parval.initcondinput(
            par.names = "Input",
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            defaults  = c(Input=0))
        par.w.expected <- streambugs:::get.inpind.parval.initcondinput(
            par.names = "w",
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            defaults  = c(w=0))
        par.initcond.input.w <- streambugs:::get.inpind.parval.initcondinput(
            par.names = c("Dini", "Input", "w"),
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            defaults  = c(Dini=0, Input=0, w=0))
        par.initcond <- streambugs:::get.single.par.inpind.parval(
            par.initcond.input.w, "Dini")
        par.input <- streambugs:::get.single.par.inpind.parval(
            par.initcond.input.w, "Input")
        par.w <- streambugs:::get.single.par.inpind.parval(
            par.initcond.input.w, "w")
        testthat::expect_equal(par.initcond, par.initcond.expected)
        testthat::expect_equal(par.input, par.input.expected)
        testthat::expect_equal(par.w, par.w.expected)
    }
})


test_that("get.inpind.parval.taxaprop", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        par.names = c("M",
                      "Ea",
                      "b",
                      "i0",
                      "EC",
                      "fbasaltax")
        par.defaults  = c(b=0.75, fbasaltax=1)
        par.taxaprop <- streambugs:::get.inpind.parval.taxaprop(
            par.names = par.names,
            y.names   = model$y.names,
            par       = model$par,
            inp       = model$inp,
            defaults  = par.defaults)
        .expect_equal.inpind.parval(
            "get_inpind_parval_taxaprop",
            par.taxaprop,
            model,
            model.name,
            par.names,
            par.defaults = par.defaults)
    }
})


test_that("get.par.proc.taxon", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        y.names = model$y.names
        par = model$par
        inp = model$inp
        # NOTE: Artificially added river width "w" env condition param to a "Miner"
        #       process to test input indices (`$inpinds`) via toy example
        kinparnames.taxon <- list(Miner    = c("kminer", "w"),
                                  Drift    = c("cdet"),
                                  Resp     = c("fresp"),
                                  Death    = c("fdeath"),
                                  Prod     = c("fprod",
                                               "fgrotax",
                                               "hdens",
                                               "KI",
                                               "KP",
                                               "KN",
                                               "w"))
        par.stoich.taxon <- streambugs:::get.par.stoich.list(
            par, names(kinparnames.taxon), streambugs:::get.par.stoich.taxon)
        # regression test
        actual.list <- streambugs:::get.par.proc.taxon(
            y.names, par, inp, kinparnames.taxon, par.stoich.taxon)
        expected.fn <- paste0("get_par_proc_taxon", "-", model.name, ".rds")
        expected.list <- readRDS(expected.fn)
        testthat::expect_equal(actual.list, expected.list)
    }
})


test_that("get.par.proc.web", {
    for (model.name in names(models)) {
        model <- models[[model.name]]
        y.names = model$y.names
        par = model$par
        inp = model$inp
        # add Pref input change to test inpinds of stoichiometrically-related taxa2
        for (inp.name in names(inp)) {
            inp[[sub("_w", "_Pref", inp.name)]] <- inp[[inp.name]]
        }
        if (!is.list(inp)) inp <- list()
        inp[paste0(y.names$habitats[1], "_Pref")] <- list(cbind(c(0, 0), c(10, 1)))
        # NOTE: Artificially added river width "w" env condition param
        #       to test input indices (`$inpinds`) via toy example
        kinparnames.web <- list(Cons     = list(taxon1 = c("fcons",
                                                           "fgrotax",
                                                           "hdens",
                                                           "Kfood",
                                                           "q",
                                                           "w"),
                                                taxa2  = c("Pref")),
                                FishPred = list(taxon1 = c("cfish",
                                                           "Kfood",
                                                           "q",
                                                           "w"),
                                                taxa2  = c("Pref")))
        par.stoich.web <- streambugs:::get.par.stoich.list(
            par, names(kinparnames.web), streambugs:::get.par.stoich.web)
        # regression test
        actual.list <- streambugs:::get.par.proc.web(
            y.names, par, inp, kinparnames.web, par.stoich.web)
        expected.fn <- paste0("get_par_proc_web", "-", model.name, ".rds")
        expected.list <- readRDS(expected.fn)
        testthat::expect_equal(actual.list, expected.list)
    }
})
