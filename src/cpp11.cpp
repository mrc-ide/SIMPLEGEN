// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// deploy_sub.cpp
cpp11::list sim_indiv_deploy(cpp11::list args, cpp11::list args_progress);
extern "C" SEXP _SIMPLEGEN_sim_indiv_deploy(SEXP args, SEXP args_progress) {
  BEGIN_CPP11
    return cpp11::as_sexp(sim_indiv_deploy(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(args), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(args_progress)));
  END_CPP11
}
// deploy_sub.cpp
cpp11::list prune_transmission_record_deploy(cpp11::list args);
extern "C" SEXP _SIMPLEGEN_prune_transmission_record_deploy(SEXP args) {
  BEGIN_CPP11
    return cpp11::as_sexp(prune_transmission_record_deploy(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(args)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_SIMPLEGEN_prune_transmission_record_deploy", (DL_FUNC) &_SIMPLEGEN_prune_transmission_record_deploy, 1},
    {"_SIMPLEGEN_sim_indiv_deploy",                 (DL_FUNC) &_SIMPLEGEN_sim_indiv_deploy,                 2},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_SIMPLEGEN(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
