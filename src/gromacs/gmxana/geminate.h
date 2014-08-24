/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef _GEMINATE_H
#define _GEMINATE_H

#include <stddef.h>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

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
                                   const real ballistic, const int nBalExp);

/* Fit to geminate recombination model.
   Returns root mean square error of fit. */
extern real fitGemRecomb(double *ct, double *time, double **ctFit,
                         const int nData, t_gemParams *params);

extern void dumpN(const real *e, const int nn, char *fn);

/* Fix NaN that might appear in the theoretical acf. */
extern void fixGemACF(double *ct, int len);

#endif
