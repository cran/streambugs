// Generated from R shell using:
//
//  > tools::package_native_routine_registration_skeleton(".")
//
// Complemented with C functions implicitly called by the deSolve package
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void streambugs_create_fA(void *, void *, void *, void *, void *);
extern void streambugs_create_fA_structure(void *);
extern void streambugs_create_input(void *, void *, void *, void *);
extern void streambugs_create_input_structure(void *);
extern void streambugs_create_parenvcondhabitat(void *, void *, void *, void *,
                                                void *, void *);
extern void streambugs_create_parenvcondhabitatgroup(void *, void *, void *,
                                                     void *, void *, void *,
                                                     void *, void *);
extern void streambugs_create_parenvcondhabitatgroups_structure(void *);
extern void streambugs_create_parenvcondreach(void *, void *, void *, void *,
                                              void *, void *);
extern void streambugs_create_parglobal(void *, void *, void *, void *, void *);
extern void streambugs_create_parglobalenvtrait(void *, void *, void *, void *,
                                                void *, void *, void *);
extern void streambugs_create_parglobalenvtraits_structure(void *);
extern void streambugs_create_parinitcond(void *, void *, void *, void *,
                                          void *, void *);
extern void streambugs_create_parinput(void *, void *, void *, void *, void *,
                                       void *);
extern void streambugs_create_partaxapropdirect(void *, void *, void *, void *,
                                                void *, void *);
extern void streambugs_create_partaxaproptrait(void *, void *, void *, void *,
                                               void *, void *, void *, void *);
extern void streambugs_create_partaxaproptraits_structure(void *);
extern void streambugs_create_processes_structure(void *);
extern void streambugs_create_proctaxon(void *, void *, void *, void *, void *,
                                        void *, void *, void *, void *, void *,
                                        void *, void *);
extern void streambugs_create_procweb(void *, void *, void *, void *, void *,
                                      void *, void *, void *);
extern void streambugs_create_procwebtaxon(void *, void *, void *, void *,
                                           void *, void *, void *, void *,
                                           void *, void *, void *, void *,
                                           void *);
extern void streambugs_delete_fA();
extern void streambugs_delete_inputs();
extern void streambugs_delete_parenvcondhabitat();
extern void streambugs_delete_parenvcondhabitatgroup();
extern void streambugs_delete_parenvcondreach();
extern void streambugs_delete_parglobal();
extern void streambugs_delete_parglobalenvtraits();
extern void streambugs_delete_parinitcond();
extern void streambugs_delete_parinput();
extern void streambugs_delete_partaxapropdirect();
extern void streambugs_delete_partaxaproptraits();
extern void streambugs_delete_processes();
extern void streambugs_init_debug(void *);

/* implicit .C calls via deSolve */
extern void streambugs_rhs_init(void *);
extern void streambugs_rhs(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"streambugs_create_fA",                                (DL_FUNC) &streambugs_create_fA,                                 5},
    {"streambugs_create_fA_structure",                      (DL_FUNC) &streambugs_create_fA_structure,                       1},
    {"streambugs_create_input",                             (DL_FUNC) &streambugs_create_input,                              4},
    {"streambugs_create_input_structure",                   (DL_FUNC) &streambugs_create_input_structure,                    1},
    {"streambugs_create_parenvcondhabitat",                 (DL_FUNC) &streambugs_create_parenvcondhabitat,                  6},
    {"streambugs_create_parenvcondhabitatgroup",            (DL_FUNC) &streambugs_create_parenvcondhabitatgroup,             8},
    {"streambugs_create_parenvcondhabitatgroups_structure", (DL_FUNC) &streambugs_create_parenvcondhabitatgroups_structure,  1},
    {"streambugs_create_parenvcondreach",                   (DL_FUNC) &streambugs_create_parenvcondreach,                    6},
    {"streambugs_create_parglobal",                         (DL_FUNC) &streambugs_create_parglobal,                          5},
    {"streambugs_create_parglobalenvtrait",                 (DL_FUNC) &streambugs_create_parglobalenvtrait,                  7},
    {"streambugs_create_parglobalenvtraits_structure",      (DL_FUNC) &streambugs_create_parglobalenvtraits_structure,       1},
    {"streambugs_create_parinitcond",                       (DL_FUNC) &streambugs_create_parinitcond,                        6},
    {"streambugs_create_parinput",                          (DL_FUNC) &streambugs_create_parinput,                           6},
    {"streambugs_create_partaxapropdirect",                 (DL_FUNC) &streambugs_create_partaxapropdirect,                  6},
    {"streambugs_create_partaxaproptrait",                  (DL_FUNC) &streambugs_create_partaxaproptrait,                   8},
    {"streambugs_create_partaxaproptraits_structure",       (DL_FUNC) &streambugs_create_partaxaproptraits_structure,        1},
    {"streambugs_create_processes_structure",               (DL_FUNC) &streambugs_create_processes_structure,                1},
    {"streambugs_create_proctaxon",                         (DL_FUNC) &streambugs_create_proctaxon,                         12},
    {"streambugs_create_procweb",                           (DL_FUNC) &streambugs_create_procweb,                            8},
    {"streambugs_create_procwebtaxon",                      (DL_FUNC) &streambugs_create_procwebtaxon,                      13},
    {"streambugs_delete_fA",                                (DL_FUNC) &streambugs_delete_fA,                                 0},
    {"streambugs_delete_inputs",                            (DL_FUNC) &streambugs_delete_inputs,                             0},
    {"streambugs_delete_parenvcondhabitat",                 (DL_FUNC) &streambugs_delete_parenvcondhabitat,                  0},
    {"streambugs_delete_parenvcondhabitatgroup",            (DL_FUNC) &streambugs_delete_parenvcondhabitatgroup,             0},
    {"streambugs_delete_parenvcondreach",                   (DL_FUNC) &streambugs_delete_parenvcondreach,                    0},
    {"streambugs_delete_parglobal",                         (DL_FUNC) &streambugs_delete_parglobal,                          0},
    {"streambugs_delete_parglobalenvtraits",                (DL_FUNC) &streambugs_delete_parglobalenvtraits,                 0},
    {"streambugs_delete_parinitcond",                       (DL_FUNC) &streambugs_delete_parinitcond,                        0},
    {"streambugs_delete_parinput",                          (DL_FUNC) &streambugs_delete_parinput,                           0},
    {"streambugs_delete_partaxapropdirect",                 (DL_FUNC) &streambugs_delete_partaxapropdirect,                  0},
    {"streambugs_delete_partaxaproptraits",                 (DL_FUNC) &streambugs_delete_partaxaproptraits,                  0},
    {"streambugs_delete_processes",                         (DL_FUNC) &streambugs_delete_processes,                          0},
    {"streambugs_init_debug",                               (DL_FUNC) &streambugs_init_debug,                                1},

    {"streambugs_rhs_init",                                 (DL_FUNC) &streambugs_rhs_init,                                  1},
    {"streambugs_rhs",                                      (DL_FUNC) &streambugs_rhs,                                       6},

    {NULL, NULL, 0}
};

void R_init_streambugs(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
