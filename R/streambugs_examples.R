################################################################
#
# streambugs 1.0
# =================
#
# ------------------
# examples of models
# ------------------
#
# creation:      07.11.2017
# modifications: 07.11.2017
#
################################################################

#' Set-up the streambugs toy model
#'
#' Set-up state variables, parameters, input, and output times of the
#' streambugs toy model. The model is ready to run with
#' \code{\link{run.streambugs}}.
#'
#' @template Streambugs_syntax
#'
#' @return List with:\describe{
#'    \item{\code{$name}}{name of the example}
#'    \item{\code{$y.names}}{list with names of state variables as returned by
#'      the \link{decode.statevarnames} function}
#'    \item{\code{$times}, \code{$par}, \code{$inp}:}{corresponding input
#'          parameters of the \code{\link{run.streambugs}} function}
#'    }
#'
#' @examples
#' model <- streambugs.example.model.toy()
#' # display values of the exponent "q" in the food limitation term; Note:
#  # non-prefixed parameter values (with "_" separator) are global values,
#  # whereas prefixed values are specific to taxa defined by the prefix
#' model$par[grepl("(.*_)?q$", names(model$par))]
#'
#' @export
streambugs.example.model.toy <- function() {

    # set-up state variable names:
    # ----------------------------

    n.Reaches         <- 3
    n.Habitats        <- 2
    # scale: reaches and habitats are independent of each-other
    n.POM             <- 2 # scale: NO (POM1 internal/autochtonous dead organic matter; POM2 unrelated external/allochtounous dead organic matter that is just mineralized or drifts away)
    n.Algae           <- 2
    n.Invertebrates   <- 5

    Reaches           <- paste("Reach",1:n.Reaches,sep="")
    Habitats          <- paste("Hab",1:n.Habitats,sep="")
    POM               <- paste("POM",1:n.POM,sep="")
    Algae             <- paste("Alga",1:n.Algae,sep="")
    Invertebrates     <- paste("Invert",1:n.Invertebrates,sep="")

    y.names <- construct.statevariables(Reaches,Habitats,POM=POM,Algae=Algae,Invertebrates=Invertebrates)
    y.names <- decode.statevarnames(y.names)


    # set-up parameters:
    # ------------------

    par <- numeric(0)

    # par["temp_vco"] = 3 + 273.15  # K
    # par["temp_col"] = 8 + 273.15  # K
    # par["temp_mod"] = 14 + 273.15  # K
    # par["temp_war"] = 22 + 273.15  # K

    # environmental parameters:

    par["w"]                              <- 10            # m
    par["Reach1_w"]                       <- 5             # m
    par["L"]                              <- 1000          # m
    par["T"]                              <- 273.15+20     # K
    par["I0"]                             <- 125           # W/m2
    par["fshade"]                         <- 0.2           # -
    par["Hab1_fshade"]                    <- 0.2           # -
    par["CP"]                             <- 0.02          # gP/m3
    par["CN"]                             <- 1             # gN/m3
    par["tau"]                            <- 0.5             # kg /(s2.m)
    par["taucrit"]                        <- 1             # kg /(s2.m)

    # initial conditions:

    par["POM_Dini"]                       <- 100           # gDM/m2
    par["Algae_Dini"]                     <- 100           # gDM/m2
    par["Invertebrates_Dini"]             <- 10            # gDM/m2

    # basal metabolism and energy content parameters:

    par["POM_EC"]                         <- 18000         # J/gDM
    par["Algae_M"]                        <- 1e-7          # g
    par["Algae_Ea"]                       <- 0.231         # eV
    par["Algae_b"]                        <- 0.93          # -
    par["Algae_i0"]                       <- 4150455552    # J/a
    par["Algae_EC"]                       <- 20000         # J/gDM
    par["Invertebrates_M"]                <- 0.001         # g
    par["Invertebrates_Ea"]               <- 0.68642       # eV
    par["Invertebrates_b"]                <- 0.695071      # -
    par["Invertebrates_i0"]               <- 1.10E+16      # J/a
    par["Invertebrates_EC"]               <- 22000         # J/gDM
    par["Algae_cdet"]                     <- 2             # s4m2/(kg2a)
    par["POM_cdet"]                       <- 3             # s4m2/(kg2a)
    par["Invertebrates_cdet"]             <- 1             # s4m2/(kg2a)

    # generate stoichiometry of mineralization, respitation, death and production:

    # POM:
    taxa.POM <- unique(y.names$y.taxa[y.names$y.groups=="POM"])
    # mineralization:
    par[paste("Miner",taxa.POM,taxa.POM,sep="_")] <- -1
    # drift:
    par[paste("Drift",taxa.POM,taxa.POM,sep="_")] <- -1



    # Algae:
    taxa.Algae <- unique(y.names$y.taxa[y.names$y.groups=="Algae"])
    # respiration:
    par[paste("Resp",taxa.Algae,taxa.Algae,sep="_")]  <- -1
    # death:
    par[paste("Death",taxa.Algae,taxa.Algae,sep="_")] <- -1
    par[paste("Death",taxa.Algae,"POM1",sep="_")]     <- 1 # scale: N/A
    # production:
    par[paste("Prod",taxa.Algae,taxa.Algae,sep="_")]  <- 1
    # drift:
    par[paste("Drift",taxa.Algae,taxa.Algae,sep="_")] <- -1


    # Invertebrates:
    taxa.Invertebrates <- unique(y.names$y.taxa[y.names$y.groups=="Invertebrates"])
    # respiration:
    par[paste("Resp",taxa.Invertebrates,taxa.Invertebrates,sep="_")]  <- -1
    # death:
    par[paste("Death",taxa.Invertebrates,taxa.Invertebrates,sep="_")] <- -1
    par[paste("Death",taxa.Invertebrates,"POM1",sep="_")]             <- 1 # scale: N/A

    # generate stoichiometry of consumption:

    par[paste("Cons",taxa.Invertebrates,"POM1",taxa.Invertebrates,sep="_")] <- 1 # scale: N/A
    par[paste("Cons",taxa.Invertebrates,"POM1","POM1",sep="_")]             <- -1 # scale: N/A
    par[paste("Cons",taxa.Invertebrates,"Alga1",taxa.Invertebrates,sep="_")] <- 1
    par[paste("Cons",taxa.Invertebrates,"Alga1","Alga1",sep="_")]             <- -1
    # scale: making more interactions, invertebrates can cons all Alga or even other
    #        invertebrates
    # scale: check if and when the steady-state is reached (preferably run until all
    #        case-studies reached steady-state) because this may affect performance
    #        if SS is reached fast scale: interactions

    # generate stoichiometry of fish predation:

    par[paste("FishPred","Fish","Invert1","Invert1",sep="_")] <- -1
    par[paste("FishPred","Fish","Invert5","Invert5",sep="_")] <- -1
    # scale: extra degradation for selected invertebrates; can make more or less of
    #        them

    # drift:
    par[paste("Drift",taxa.Invertebrates,taxa.Invertebrates,sep="_")] <- -1


    # process parameters:

    par["POM_kminer"]                     <- 20            # 1/a

    par["Algae_hdens"]                    <- 1000          # gDM/m2
    par["Algae_KI"]                       <- 62.1          # W/m2
    par["Algae_KP"]                       <- 0.01          # gP/m3
    par["Algae_KN"]                       <- 0.1           # gN/m3
    par["Algae_fresp"]                    <- 1             # -
    par["Algae_fdeath"]                   <- 0.7           # -
    par["Algae_fprod"]                    <- 5             # -

    par["Invertebrates_hdens"]            <- 10            # gDM/m2
    par["Invertebrates_fresp"]            <- 2.5           # -
    par["Invertebrates_fdeath"]           <- 0.7           # -
    par["Invertebrates_fcons"]            <- 5             # -
    par["Invertebrates_Kfood"]            <- 1             # gDM/m2
    par["Invertebrates_q"]                <- 2
    par["Invertebrates_Pref"]             <- 1             # -

    par["DFish"]                          <- 100           # kg/ha
    par["Reach2_DFish"]                   <- 200           # kg/ha
    # scale: rm all reach/habitat-dependent parameters and make reaches and habitat
    #        homogeneous (here: not _d param is a value for all)
    par["Fish_cfish"]                     <- 10            # gDM/kg/d
    par["Fish_Kfood"]                     <- 1             # gDM/m2
    par["Fish_q"]                         <- 1             # -
    par["Fish_Pref"]                      <- 1             # -

    # set-up inputs:
    # --------------

    inp <- list(Reach3_w=matrix(c(0:1,10,2),ncol=2,byrow=FALSE))
    # scale: use no inputs; here: width of Reach3 (river) decreases at t=1 from 10
    # to 2 (2nd row)
    #inp <- NA

    # set-up output times:
    # --------------------

    tout <- 0:200/100

    return(list(
        name    = "Toy example",
        y.names = y.names,
        times   = tout,
        par     = par,
        inp     = inp
    ))
}
