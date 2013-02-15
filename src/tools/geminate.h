#ifndef _GEMINATE_H
#define _GEMINATE_H

enum {
    gemNULL, gemNONE, gemDD, gemAD, gemAA, gemA4, gemNR
};
static const char *gemType[] = {NULL, "none", "dd", "ad", "aa", "a4", NULL};

/* The first few sections of this file contain functions that were adopted,
 * and to some extent modified, by Erik Marklund (erikm[aT]xray.bmc.uu.se,
 * http://folding.bmc.uu.se) from code written by Omer Markovitch (email, url).
 * This is also the case with the function eq10v2() in geminate.c.
 *
 * The parts menetioned in the previous paragraph were contributed under a BSD license.
 */

/* This first part is derived from complex.c which I recieved from Omer Markowitch.
 * - Erik Marklund
 *
 * ------------- from complex.c ------------- */

#include <math.h>
/* definition of PI */
#ifndef PI
#define PI (acos(-1.0))
#endif

/* definition of the type complex */
typedef struct
{
    double r, i;
} gem_complex;


/* ------------ end of complex.c ------------ */

/* This next part was derived from cerror.c and rerror.c,
 * also received from Omer Markovitch.
 * ------------- from [cr]error.c ------------- */

#ifndef sPI
#define sPI (sqrt(PI))
#endif

/* ------------ end of [cr]error.c ------------ */

/* ///////////////// REVERSIBLE GEMINATE RECOMBINATION ///////////////////
 * Here follow routines and structs for reversible geminate recombination.
 */

typedef struct {
    size_t  n;
    double *y;
    double  tDelta;
    int     nexp;
} balData;


typedef struct {
    /* Used in the rewritten version of Omer's gem. recomb. analysis */
    double ka, kd;               /* -- || -- results[]  */
    double sigma;                /* -- || -- sigma      */
    double D;                    /* The diffusion coefficient */
    double kD;                   /* Replaces kD in analyse_corr_gem3d() */

    /* The following members are for calcsquare() and takeAwayBallistic() */
    double tDelta;            /* Time between frames */
    /* double logAfterTime;        /\* Time after which we do the lsq calculations on a logarithmic timescale. *\/ */
    int    nFitPoints;        /* Number of points to take from the ACF for fitting */
    double begFit;            /* Fit from this time (ps) */
    double endFit;            /* Fit up to this time (ps) */
/*   double logDelta; */
/*   double logPF; */
/* To get an equal number of points in the lin and log regime,
 * we'll use logDelta to calculate where to sample the ACF.
 * if i and j are indices in the linear and log regime, then:
 *   j = Cexp(A(i+nLin)),
 * where C = (nLin**2 / len) and A = log(len/nLin) / nLin.
 * This expands to j = (nLin**2 / len) * exp((i+nLin) * log(len/nLin) / nLin).
 * We'll save part of our brains and some computing time if we pre-calculate
 *  1) logDelta = log(len/nLin) / nLin
 *  2) logPF    = nLin**2 / len
 * and get j = logPF * exp((i+nLin) * logDelta). */

    /* I will redo this for a fit done entirely in log-log.
     *  j' = j+1
     *  nFitPoints' = nFitPoints-1
     *
     *  j' = Cexp(Ai)
     *  (j'= 1 <=> i=0)
     *     => C=1
     *  (j'=len <=> i=nFitPoints')
     *     => A=log(len)/nFitPoints'
     *     => j = exp(i*log(len)/(nFitPoints-1)) -1
     **/
/* #define GETLOGINDEX(i,params) (params)->logPF * exp(((i)+(params)->nLin) * (params)->logDelta)
 */
    double   logQuota;
    int      nLin;          /* Number of timepoints in the linear regime */
    int      len;           /* Length of time and ct arrays */
    int      nExpFit;       /* Number of exponentials to fit */
    real     ballistic;     /* Time before which the ballistic term should be fitted */
    gmx_bool bDt;           /* TRUE =>  use time derivative at time 0
                             *          to find fastest component.
                             * FALSE => use coefficient in exponenetial
                             *          to find fastest component. */
} t_gemParams;


typedef struct {
    size_t       n;        /* Number of data points (lin-log) */
    double      *y;        /* The datapoints */
    double      *ctTheory; /* The theoretical ct which will be built by gemFunc_f. */
    double      *LinLog;
    double      *time;
    double       ka;
    double       kd;
    double       tDelta; /* time difference between subsequent datapoints */
    size_t       nData;  /* real size of the data */
    int         *logtime;
    double      *doubleLogTime;
    t_gemParams *params;
} gemFitData;

extern void takeAwayBallistic(double *ct, double *t,
                              int len, real tMax,
                              int nexp, gmx_bool bDerivative);


extern t_gemParams *init_gemParams(const double sigma, const double D,
                                   const real *t, const int len, const int nFitPoints,
                                   const real begFit, const real endFit,
                                   const real ballistic, const int nBalExp, const gmx_bool bDt);

/* Fit to geminate recombination model.
   Returns root mean square error of fit. */
extern real fitGemRecomb(double *ct, double *time, double **ctFit,
                         const int nData, t_gemParams *params);

extern void dumpN(const real *e, const int nn, char *fn);

/* Fix NaN that might appear in the theoretical acf. */
extern void fixGemACF(double *ct, int len);

#endif
