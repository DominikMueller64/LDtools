#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP LDtools_center_matrix(SEXP, SEXP);
extern SEXP LDtools_comb_adj(SEXP);
extern SEXP LDtools_comb_all(SEXP);
extern SEXP LDtools_comb_all_sets(SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_flank_sets(SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_nearest_k_sets(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_sliding(SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_sliding_sets(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_wind(SEXP, SEXP, SEXP);
extern SEXP LDtools_comb_wind_sets(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_find_closest(SEXP, SEXP);
extern SEXP LDtools_get_counts(SEXP, SEXP);
extern SEXP LDtools_index_geq(SEXP, SEXP);
extern SEXP LDtools_index_greater(SEXP, SEXP);
extern SEXP LDtools_LD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_LD_mult(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_LD_mult_r(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_LD_mult_r_dev(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_optimum_pxy(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP LDtools_sliding_window(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"LDtools_center_matrix",       (DL_FUNC) &LDtools_center_matrix,       2},
    {"LDtools_comb_adj",            (DL_FUNC) &LDtools_comb_adj,            1},
    {"LDtools_comb_all",            (DL_FUNC) &LDtools_comb_all,            1},
    {"LDtools_comb_all_sets",       (DL_FUNC) &LDtools_comb_all_sets,       4},
    {"LDtools_comb_flank_sets",     (DL_FUNC) &LDtools_comb_flank_sets,     4},
    {"LDtools_comb_nearest_k_sets", (DL_FUNC) &LDtools_comb_nearest_k_sets, 5},
    {"LDtools_comb_sliding",        (DL_FUNC) &LDtools_comb_sliding,        4},
    {"LDtools_comb_sliding_sets",   (DL_FUNC) &LDtools_comb_sliding_sets,   7},
    {"LDtools_comb_wind",           (DL_FUNC) &LDtools_comb_wind,           3},
    {"LDtools_comb_wind_sets",      (DL_FUNC) &LDtools_comb_wind_sets,      6},
    {"LDtools_find_closest",        (DL_FUNC) &LDtools_find_closest,        2},
    {"LDtools_get_counts",          (DL_FUNC) &LDtools_get_counts,          2},
    {"LDtools_index_geq",           (DL_FUNC) &LDtools_index_geq,           2},
    {"LDtools_index_greater",       (DL_FUNC) &LDtools_index_greater,       2},
    {"LDtools_LD",                  (DL_FUNC) &LDtools_LD,                  7},
    {"LDtools_LD_mult",             (DL_FUNC) &LDtools_LD_mult,             5},
    {"LDtools_LD_mult_r",           (DL_FUNC) &LDtools_LD_mult_r,           5},
    {"LDtools_LD_mult_r_dev",       (DL_FUNC) &LDtools_LD_mult_r_dev,       6},
    {"LDtools_optimum_pxy",         (DL_FUNC) &LDtools_optimum_pxy,         5},
    {"LDtools_sliding_window",      (DL_FUNC) &LDtools_sliding_window,      6},
    {NULL, NULL, 0}
};

void R_init_LDtools(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
