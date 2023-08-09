#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(kdenestmlcv)(void *, void *, void *, void *, void *);
extern void F77_NAME(kdenestmlcvb)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srhointegrand)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srhointegrand2)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srhointegrandv)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(srhosum)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ssbiv)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ssbiv2)(void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ssbivb)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ssuni)(void *, void *, void *, void *, void *);
extern void F77_NAME(ssuni2)(void *, void *, void *, void *, void *);
extern void F77_NAME(ssunib)(void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(surrogateacf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"kdenestmlcv",   (DL_FUNC) &F77_NAME(kdenestmlcv),    5},
    {"kdenestmlcvb",  (DL_FUNC) &F77_NAME(kdenestmlcvb),   6},
    {"srhointegrand", (DL_FUNC) &F77_NAME(srhointegrand),  9},
    {"srhointegrand2",(DL_FUNC) &F77_NAME(srhointegrand2), 8},
    {"srhointegrandv",(DL_FUNC) &F77_NAME(srhointegrandv),10},
    {"srhosum",       (DL_FUNC) &F77_NAME(srhosum),        8},
    {"ssbiv",         (DL_FUNC) &F77_NAME(ssbiv),          6},
    {"ssbiv2",        (DL_FUNC) &F77_NAME(ssbiv2),         6},
    {"ssbivb",        (DL_FUNC) &F77_NAME(ssbivb),         9},
    {"ssuni",         (DL_FUNC) &F77_NAME(ssuni),          5},
    {"ssuni2",        (DL_FUNC) &F77_NAME(ssuni2),         5},
    {"ssunib",        (DL_FUNC) &F77_NAME(ssunib),         8},
    {"surrogateacf",  (DL_FUNC) &F77_NAME(surrogateacf),  11},
    {NULL, NULL, 0}
};

void R_init_tseriesEntropy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
