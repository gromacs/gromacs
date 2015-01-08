/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal
 * \file
 * \brief
 * Implements routine for fitting a data set to a curve
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "expfit.h"

#include <math.h>
#include <string.h>

#include <algorithm>

#include "external/lmfit/lmcurve.h"

#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Number of parameters for each fitting function */
static const int nfp_ffn[effnNR] = { 0, 1, 2, 3, 5, 7, 9, 2, 4, 3, 6 };

/* +2 becuase parse_common_args wants leading and concluding NULL.
 * We only allow exponential functions as choices on the command line,
 * hence there are many more NULL field (which have to be at the end of
 * the array).
 */
const char      *s_ffn[effnNR+2] = {
    NULL, "none", "exp", "aexp", "exp_exp",
    "exp5", "exp7", "exp9",
    NULL, NULL, NULL, NULL, NULL
};

/*! \brief Long description for each fitting function type */
static const char *longs_ffn[effnNR] = {
    "no fit",
    "y = exp(-x/|a0|)",
    "y = a1 exp(-x/|a0|)",
    "y = a1 exp(-x/|a0|) + (1-a1) exp(-x/(|a2|)), a2 > a0 > 0",
    "y = a0 exp(-x/|a1|) +  a2 exp(-x/|a3|) + a4, a3 >= a1",
    "y = a0 exp(-x/|a1|) +  a2 exp(-x/|a3|) + a4 exp(-x/|a5|) + a6, a5 >= a3 >= a1",
    "y = a0 exp(-x/|a1|) +  a2 exp(-x/|a3|) + a4 exp(-x/|a5|) + a6 exp(-x/|a7|) + a8, a7 >= a5 >= a3 >= a1",
    "y = exp(-v) (cosh(wv) + 1/w sinh(wv)), v = x/(2 a0), w = sqrt(1 - a1)",
    "y = 1/2*(a0+a1) - 1/2*(a0-a1)*erf( (x-a2) /a3^2)",
    "y = a1 *2*a0*((exp(-x/a0)-1)*(a0/x)+1)+(1-a1)*2*a2*((exp(-x/a2)-1)*(a2/x)+1)",
    "y = (1-a0)*cos(x*a1)*exp(-(x/a2)^a3) + a0*exp(-(x/a4)^a5)"
};

int effnNparams(int effn)
{
    if ((0 <= effn) && (effn < effnNR))
    {
        return nfp_ffn[effn];
    }
    else
    {
        return -1;
    }
}

const char *effnDescription(int effn)
{
    if ((0 <= effn) && (effn < effnNR))
    {
        return longs_ffn[effn];
    }
    else
    {
        return NULL;
    }
}

int sffn2effn(const char **sffn)
{
    int eFitFn, i;

    eFitFn = 0;
    for (i = 0; i < effnNR; i++)
    {
        if (sffn[i+1] && strcmp(sffn[0], sffn[i+1]) == 0)
        {
            eFitFn = i;
        }
    }

    return eFitFn;
}

/*! \brief Compute exponential function A exp(-x/tau) */
static double myexp(double x, double A, double tau)
{
    if ((A == 0) || (tau == 0))
    {
        return 0;
    }
    return A*exp(-x/tau);
}

/*! \brief Compute y=(a0+a1)/2-(a0-a1)/2*erf((x-a2)/a3^2) */
static double lmc_erffit (double x, const double *a)
{
    double erfarg;
    double myerf;

    if (a[3] != 0)
    {
        erfarg = (x-a[2])/(a[3]*a[3]);
        myerf  = gmx_erfd(erfarg);
    }
    else
    {
        /* If a[3] == 0, a[3]^2 = 0 and the erfarg becomes +/- infinity */
        if (x < a[2])
        {
            myerf = -1;
        }
        else
        {
            myerf = 1;
        }
    }
    return 0.5*((a[0]+a[1]) - (a[0]-a[1])*myerf);
}

/*! \brief Exponent function that prevents overflow */
static double safe_exp(double x)
{
    double exp_max = 200;
    double exp_min = -exp_max;
    if (x <= exp_min)
    {
        return exp(exp_min);
    }
    else if (x >= exp_max)
    {
        return exp(exp_max);
    }
    else
    {
        return exp(x);
    }
}

/*! \brief Exponent minus 1 function that prevents overflow */
static double safe_expm1(double x)
{
    double exp_max = 200;
    double exp_min = -exp_max;
    if (x <= exp_min)
    {
        return -1;
    }
    else if (x >= exp_max)
    {
        return exp(exp_max);
    }
    else
    {
        return gmx_expm1(x);
    }
}

/*! \brief Compute y = exp(-x/|a0|) */
static double lmc_exp_one_parm(double x, const double *a)
{
    return safe_exp(-x/fabs(a[0]));
}

/*! \brief Compute y = a1 exp(-x/|a0|) */
static double lmc_exp_two_parm(double x, const double *a)
{
    return a[1]*safe_exp(-x/fabs(a[0]));
}

/*! \brief Compute y = a1 exp(-x/|a0|) + (1-a1) exp(-x/|a2|) */
static double lmc_exp_exp(double x, const double *a)
{
    double e1, e2;

    e1      = safe_exp(-x/fabs(a[0]));
    e2      = safe_exp(-x/(fabs(a[0])+fabs(a[2])));
    return a[1]*e1 + (1-a[1])*e2;
}

/*! \brief Compute y = a0 exp(-x/|a1|) + a2 exp(-x/(|a1|+|a3|)) + a4 */
static double lmc_exp_5_parm(double x, const double *a)
{
    double e1, e2;

    e1      = safe_exp(-x/fabs(a[1]));
    e2      = safe_exp(-x/(fabs(a[1])+fabs(a[3])));
    return a[0]*e1 + a[2]*e2 + a[4];
}

/*! \brief Compute 7 parameter exponential function value.
 *
 * Compute y = a0 exp(-x/|a1|) + a2 exp(-x/(|a1|+|a3|)) +
 * a4 exp(-x/(|a1|+|a3|+|a5|)) + a6
 */
static double lmc_exp_7_parm(double x, const double *a)
{
    double e1, e2, e3;
    double fa1, fa3, fa5;

    fa1 = fabs(a[1]);
    fa3 = fa1 + fabs(a[3]);
    fa5 = fa3 + fabs(a[5]);
    e1  = safe_exp(-x/fa1);
    e2  = safe_exp(-x/fa3);
    e3  = safe_exp(-x/fa5);
    return a[0]*e1 + a[2]*e2 + a[4]*e3 + a[6];
}

/*! \brief Compute 9 parameter exponential function value.
 *
 * Compute y = a0 exp(-x/|a1|) + a2 exp(-x/(|a1|+|a3|)) +
 * a4 exp(-x/(|a1|+|a3|+|a5|)) + a6 exp(-x/(|a1|+|a3|+|a5|+|a7|)) + a8
 */
static double lmc_exp_9_parm(double x, const double *a)
{
    double e1, e2, e3, e4;
    double fa1, fa3, fa5, fa7;

    fa1 = fabs(a[1]);
    fa3 = fa1 + fabs(a[3]);
    fa5 = fa3 + fabs(a[5]);
    fa7 = fa5 + fabs(a[7]);

    e1      = safe_exp(-x/fa1);
    e2      = safe_exp(-x/fa3);
    e3      = safe_exp(-x/fa5);
    e4      = safe_exp(-x/fa7);
    return a[0]*e1 + a[2]*e2 + a[4]*e3 + a[6]*e4 + a[8];
}

/*! \brief Compute y = (1-a0)*exp(-(x/|a2|)^|a3|)*cos(x*|a1|) + a0*exp(-(x/|a4|)^|a5|) */
static double lmc_pres_6_parm(double x, const double *a)
{
    double term1, term2, term3;
    double pow_max = 10;

    term3  = 0;
    if ((a[4] != 0) && (a[0] != 0))
    {
        double power = std::min(fabs(a[5]), pow_max);
        term3        = a[0] * safe_exp(-pow((x/fabs(a[4])), power));
    }

    term1  = 1-a[0];
    term2  = 0;
    if ((term1 != 0) && (a[2] != 0))
    {
        double power = std::min(fabs(a[3]), pow_max);
        term2        = safe_exp(-pow((x/fabs(a[2])), power)) * cos(x*fabs(a[1]));
    }

    return term1*term2 + term3;
}

/*! \brief Compute vac function */
static double lmc_vac_2_parm(double x, const double *a)
{
    /* Fit to function
     *
     * y = 1/2 (1 - 1/w) exp(-(1+w)v) + 1/2 (1 + 1/w) exp(-(1-w)v)
     *
     *   = exp(-v) (cosh(wv) + 1/w sinh(wv))
     *
     *    v = x/(2 a0)
     *    w = sqrt(1 - a1)
     *
     *    For tranverse current autocorrelation functions:
     *       a0 = tau
     *       a1 = 4 tau (eta/rho) k^2
     *
     */

    double y, v, det, omega, wv, em, ec, es;
    double wv_max = 100;

    v   = x/(2*fabs(a[0]));
    det = 1 - a[1];
    em  = safe_exp(-v);
    if (det != 0)
    {
        omega = sqrt(fabs(det));
        wv    = std::min(omega*v, wv_max);

        if (det > 0)
        {
            ec = em*0.5*(safe_exp(wv)+safe_exp(-wv));
            es = em*0.5*(safe_exp(wv)-safe_exp(-wv))/omega;
        }
        else
        {
            ec = em*cos(wv);
            es = em*sin(wv)/omega;
        }
        y      = ec + es;
    }
    else
    {
        y      = (1+v)*em;
    }
    return y;
}

/*! \brief Compute error estimate */
static double lmc_errest_3_parm(double x, const double *a)
{
    double e1, e2, v1, v2;
    double fa0 = fabs(a[0]);
    double fa1;
    double fa2 = fa0+fabs(a[2]);

    if (a[0] != 0)
    {
        e1 = safe_expm1(-x/fa0);
    }
    else
    {
        e1 = 0;
    }
    if (a[2] != 0)
    {
        e2 = safe_expm1(-x/fa2);
    }
    else
    {
        e2 = 0;
    }

    if (x > 0)
    {
        v1      = 2*fa0*(e1*fa0/x + 1);
        v2      = 2*fa2*(e2*fa2/x + 1);
        /* We need 0 <= a1 <= 1 */
        fa1     = std::min(1.0, std::max(0.0, a[1]));

        return fa1*v1 + (1-fa1)*v2;
    }
    else
    {
        return 0;
    }
}

/*! \brief function type for passing to fitting routine */
typedef double (*t_lmcurve)(double x, const double *a);

/*! \brief array of fitting functions corresponding to the pre-defined types */
t_lmcurve lmcurves[effnNR+1] = {
    lmc_exp_one_parm, lmc_exp_one_parm, lmc_exp_two_parm,
    lmc_exp_exp,      lmc_exp_5_parm,   lmc_exp_7_parm,
    lmc_exp_9_parm,
    lmc_vac_2_parm,   lmc_erffit,  lmc_errest_3_parm, lmc_pres_6_parm
};

double fit_function(const int eFitFn, const double parm[], const double x)
{
    if ((eFitFn < 0) || (eFitFn >= effnNR))
    {
        fprintf(stderr, "fitfn = %d, should be in the range 0..%d\n",
                eFitFn, effnNR-1);
        return 0.0;
    }
    return lmcurves[eFitFn](x, parm);
}

/*! \brief lmfit_exp supports fitting of different functions
 *
 * This routine calls the Levenberg-Marquardt non-linear fitting
 * routine for fitting a data set with errors to a target function.
 * Fitting routines included in gromacs in src/external/lmfit.
 */
static gmx_bool lmfit_exp(int          nfit,
                          const double x[],
                          const double y[],
                          const double dy[],
                          double       parm[],
                          gmx_bool     bVerbose,
                          int          eFitFn,
                          int          nfix)
{
    double             chisq, ochisq;
    gmx_bool           bCont;
    int                j;
    int                maxiter = 100;
    lm_control_struct  control;
    lm_status_struct  *status;
    int                nparam = effnNparams(eFitFn);
    int                p2;
    gmx_bool           bSkipLast;

    if ((eFitFn < 0) || (eFitFn >= effnNR))
    {
        fprintf(stderr, "fitfn = %d, should be in the range 0..%d\n",
                eFitFn, effnNR-1);
        return FALSE;
    }
    /* Using default control structure for double precision fitting that
     * comes with the lmfit package (i.e. from the include file).
     */
    control            = lm_control_double;
    control.verbosity  = (bVerbose ? 1 : 0);
    control.n_maxpri   = 0;
    control.m_maxpri   = 0;

    snew(status, 1);
    /* Initial params */
    chisq  = 1e12;
    j      = 0;
    if (bVerbose)
    {
        printf("%4s  %10s  Parameters\n", "Step", "chi^2");
    }
    /* Check whether we have to skip some params */
    if (nfix > 0)
    {
        do
        {
            p2        = 1 << (nparam-1);
            bSkipLast = ((p2 & nfix) == p2);
            if (bSkipLast)
            {
                nparam--;
                nfix -= p2;
            }
        }
        while ((nparam > 0) && (bSkipLast));
        if (bVerbose)
        {
            printf("Using %d out of %d parameters\n", nparam, effnNparams(eFitFn));
        }
    }
    do
    {
        ochisq = chisq;
        lmcurve(nparam, parm, nfit, x, y, dy,
                lmcurves[eFitFn], &control, status);
        chisq = sqr(status->fnorm);
        if (bVerbose)
        {
            printf("status: fnorm = %g, nfev = %d, userbreak = %d\noutcome = %s\n",
                   status->fnorm, status->nfev, status->userbreak,
                   lm_infmsg[status->outcome]);
        }
        if (bVerbose)
        {
            int mmm;
            printf("%4d  %8g", j, chisq);
            for (mmm = 0; (mmm < effnNparams(eFitFn)); mmm++)
            {
                printf("  %8g", parm[mmm]);
            }
            printf("\n");
        }
        j++;
        bCont = (fabs(ochisq - chisq) > fabs(control.ftol*chisq));
    }
    while (bCont && (j < maxiter));

    sfree(status);

    return TRUE;
}

/*! \brief Ensure the fitting parameters are well-behaved.
 *
 * In order to make sure that fitting behaves according to the
 * constraint that time constants are positive and increasing
 * we enforce this here by setting all time constants to their
 * absolute value and by adding e.g. |a_0| to |a_2|. This is
 * done in conjunction with the fitting functions themselves.
 * When there are multiple time constants we make sure that
 * the successive time constants are at least double the
 * previous ones and during fitting we enforce the they remain larger.
 * It may very well help the convergence of the fitting routine.
 */
static void initiate_fit_params(int    eFitFn,
                                double params[])
{
    int i, nparm;

    nparm = effnNparams(eFitFn);

    switch (eFitFn)
    {
        case effnVAC:
            GMX_ASSERT(params[0] >= 0, "parameters should be >= 0");
            break;
        case effnEXP1:
        case effnEXP2:
        case effnEXPEXP:
            GMX_ASSERT(params[0] >= 0, "parameters should be >= 0");
            if (nparm > 2)
            {
                GMX_ASSERT(params[2] >= 0, "parameters should be >= 0");
                /* In order to maintain params[2] >= params[0] in the final
                 * result, we fit the difference between the two, that is
                 * params[2]-params[0] and in the add add in params[0]
                 * again.
                 */
                params[2] = std::max(fabs(params[2])-params[0], params[0]);
            }
            break;
        case effnEXP5:
        case effnEXP7:
        case effnEXP9:
            GMX_ASSERT(params[1] >= 0, "parameters should be >= 0");
            params[1] = fabs(params[1]);
            if (nparm > 3)
            {
                GMX_ASSERT(params[3] >= 0, "parameters should be >= 0");
                /* See comment under effnEXPEXP */
                params[3] = std::max(fabs(params[3])-params[1], params[1]);
                if (nparm > 5)
                {
                    GMX_ASSERT(params[5] >= 0, "parameters should be >= 0");
                    /* See comment under effnEXPEXP */
                    params[5] = std::max(fabs(params[5])-params[3], params[3]);
                    if (nparm > 7)
                    {
                        GMX_ASSERT(params[7] >= 0, "parameters should be >= 0");
                        /* See comment under effnEXPEXP */
                        params[7] = std::max(fabs(params[7])-params[5], params[5]);
                    }
                }
            }
            break;
        case effnERREST:
            GMX_ASSERT(params[0] >= 0, "parameters should be >= 0");
            GMX_ASSERT(params[1] >= 0 && params[1] <= 1, "parameter 1 should in 0 .. 1");
            GMX_ASSERT(params[2] >= 0, "parameters should be >= 0");
            /* See comment under effnEXPEXP */
            params[2] = fabs(params[2])-params[0];
            break;
        case effnPRES:
            for (i = 1; (i < nparm); i++)
            {
                GMX_ASSERT(params[i] >= 0, "parameters should be >= 0");
            }
            break;
        default:
            break;
    }
}

/*! \brief Process the fitting parameters to get output parameters.
 *
 * See comment at the previous function.
 */
static void extract_fit_params(int    eFitFn,
                               double params[])
{
    int i, nparm;

    nparm = effnNparams(eFitFn);

    switch (eFitFn)
    {
        case effnVAC:
            params[0] = fabs(params[0]);
            break;
        case effnEXP1:
        case effnEXP2:
        case effnEXPEXP:
            params[0] = fabs(params[0]);
            if (nparm > 2)
            {
                /* Back conversion of parameters from the fitted difference
                 * to the absolute value.
                 */
                params[2] = fabs(params[2])+params[0];
            }
            break;
        case effnEXP5:
        case effnEXP7:
        case effnEXP9:
            params[1] = fabs(params[1]);
            if (nparm > 3)
            {
                /* See comment under effnEXPEXP */
                params[3] = fabs(params[3])+params[1];
                if (nparm > 5)
                {
                    /* See comment under effnEXPEXP */
                    params[5] = fabs(params[5])+params[3];
                    if (nparm > 7)
                    {
                        /* See comment under effnEXPEXP */
                        params[7] = fabs(params[7])+params[5];
                    }
                }
            }
            break;
        case effnERREST:
            params[0] = fabs(params[0]);
            if (params[1] < 0)
            {
                params[1] = 0;
            }
            else if (params[1] > 1)
            {
                params[1] = 1;
            }
            /* See comment under effnEXPEXP */
            params[2] = params[0]+fabs(params[2]);
            break;
        case effnPRES:
            for (i = 1; (i < nparm); i++)
            {
                params[i] = fabs(params[i]);
            }
            break;
        default:
            break;
    }
}

/*! \brief Print chi-squared value and the parameters */
static void print_chi2_params(FILE        *fp,
                              const int    eFitFn,
                              const double fitparms[],
                              const char  *label,
                              const int    nfitpnts,
                              const double x[],
                              const double y[])
{
    int    i;
    double chi2 = 0;

    for (i = 0; (i < nfitpnts); i++)
    {
        double yfit = lmcurves[eFitFn](x[i], fitparms);
        chi2 += sqr(y[i] - yfit);
    }
    fprintf(fp, "There are %d data points, %d parameters, %s chi2 = %g\nparams:",
            nfitpnts, effnNparams(eFitFn), label, chi2);
    for (i = 0; (i < effnNparams(eFitFn)); i++)
    {
        fprintf(fp, "  %10g", fitparms[i]);
    }
    fprintf(fp, "\n");
}

/*! \brief See description in header file. */
real do_lmfit(int ndata, real c1[], real sig[], real dt, real x0[],
              real begintimefit, real endtimefit, const output_env_t oenv,
              gmx_bool bVerbose, int eFitFn, double fitparms[], int fix,
              const char *fn_fitted)
{
    FILE    *fp;
    int      i, j, nfitpnts;
    double   integral, ttt;
    double  *x, *y, *dy;

    if (0 != fix)
    {
        fprintf(stderr, "Using fixed parameters in curve fitting is temporarily not working.\n");
    }
    if (debug)
    {
        fprintf(debug, "There are %d points to fit %d vars!\n", ndata, effnNparams(eFitFn));
        fprintf(debug, "Fit to function %d from %g through %g, dt=%g\n",
                eFitFn, begintimefit, endtimefit, dt);
    }

    snew(x, ndata);
    snew(y, ndata);
    snew(dy, ndata);
    j    = 0;
    for (i = 0; i < ndata; i++)
    {
        ttt = x0 ? x0[i] : dt*i;
        if (ttt >= begintimefit && ttt <= endtimefit)
        {
            x[j]  = ttt;
            y[j]  = c1[i];
            if (NULL == sig)
            {
                // No weighting if all values are divided by 1.
                dy[j] = 1;
            }
            else
            {
                dy[j] = std::max(1.0e-7, (double)sig[i]);
            }
            if (debug)
            {
                fprintf(debug, "j= %d, i= %d, x= %g, y= %g, dy=%g, ttt=%g\n",
                        j, i, x[j], y[j], dy[j], ttt);
            }
            j++;
        }
    }
    nfitpnts = j;
    integral = 0;
    if (nfitpnts < effnNparams(eFitFn))
    {
        fprintf(stderr, "Not enough (%d) data points for fitting, dt = %g!\n",
                nfitpnts, dt);
    }
    else
    {
        gmx_bool bSuccess;

        if (bVerbose)
        {
            print_chi2_params(stdout, eFitFn, fitparms, "initial", nfitpnts, x, y);
        }
        initiate_fit_params(eFitFn, fitparms);

        bSuccess = lmfit_exp(nfitpnts, x, y, dy, fitparms, bVerbose, eFitFn, fix);
        extract_fit_params(eFitFn, fitparms);

        if (!bSuccess)
        {
            fprintf(stderr, "Fit failed!\n");
        }
        else
        {
            if (bVerbose)
            {
                print_chi2_params(stdout, eFitFn, fitparms, "final", nfitpnts, x, y);
            }
            switch (eFitFn)
            {
                case effnEXP1:
                    integral = fitparms[0]*myexp(begintimefit, 1,  fitparms[0]);
                    break;
                case effnEXP2:
                    integral = fitparms[0]*myexp(begintimefit, fitparms[1],  fitparms[0]);
                    break;
                case effnEXPEXP:
                    integral = (fitparms[0]*myexp(begintimefit, fitparms[1],  fitparms[0]) +
                                fitparms[2]*myexp(begintimefit, 1-fitparms[1], fitparms[2]));
                    break;
                case effnEXP5:
                case effnEXP7:
                case effnEXP9:
                    integral = 0;
                    for (i = 0; (i < (effnNparams(eFitFn)-1)/2); i++)
                    {
                        integral += fitparms[2*i]*myexp(begintimefit, fitparms[2*i+1], fitparms[2*i]);
                    }
                    break;
                default:
                    /* Do numerical integration */
                    integral = 0;
                    for (i = 0; (i < nfitpnts-1); i++)
                    {
                        double y0 = lmcurves[eFitFn](x[i], fitparms);
                        double y1 = lmcurves[eFitFn](x[i+1], fitparms);
                        integral += (x[i+1]-x[i])*(y1+y0)*0.5;
                    }
                    break;
            }

            if (bVerbose)
            {
                printf("FIT: Integral of fitted function: %g\n", integral);
                if ((effnEXP5 == eFitFn) || (effnEXP7 == eFitFn) || (effnEXP9 == eFitFn))
                {
                    printf("FIT: Note that the constant term is not taken into account when computing integral.\n");
                }
            }
            /* Generate debug output */
            if (NULL != fn_fitted)
            {
                fp = xvgropen(fn_fitted, "Data + Fit", "Time (ps)",
                              "Data (t)", oenv);
                for (i = 0; (i < effnNparams(eFitFn)); i++)
                {
                    fprintf(fp, "# fitparms[%d] = %g\n", i, fitparms[i]);
                }
                for (j = 0; (j < nfitpnts); j++)
                {
                    real ttt = x0 ? x0[i] : dt*j;
                    fprintf(fp, "%10.5e  %10.5e  %10.5e\n",
                            x[j], y[j], lmcurves[eFitFn](ttt, fitparms));
                }
                xvgrclose(fp);
            }
        }
    }

    sfree(x);
    sfree(y);
    sfree(dy);

    return integral;
}

/*! See description in header file. */
real fit_acf(int ncorr, int fitfn, const output_env_t oenv, gmx_bool bVerbose,
             real tbeginfit, real tendfit, real dt, real c1[], real *fit)
{
    double      fitparm[3];
    double      tStart, tail_corr, sum, sumtot = 0, c_start, ct_estimate;
    real       *sig;
    int         i, j, jmax, nf_int;
    gmx_bool    bPrint;

    bPrint = bVerbose || bDebugMode();

    if (bPrint)
    {
        printf("COR:\n");
    }

    if (tendfit <= 0)
    {
        tendfit = ncorr*dt;
    }
    nf_int = std::min(ncorr, (int)(tendfit/dt));
    sum    = print_and_integrate(debug, nf_int, dt, c1, NULL, 1);

    /* Estimate the correlation time for better fitting */
    ct_estimate = 0.5*c1[0];
    for (i = 1; (i < ncorr) && (c1[i] > 0); i++)
    {
        ct_estimate += c1[i];
    }
    ct_estimate *= dt/c1[0];

    if (bPrint)
    {
        printf("COR: Correlation time (plain integral from %6.3f to %6.3f ps) = %8.5f ps\n",
               0.0, dt*nf_int, sum);
        printf("COR: Relaxation times are computed as fit to an exponential:\n");
        printf("COR:   %s\n", effnDescription(fitfn));
        printf("COR: Fit to correlation function from %6.3f ps to %6.3f ps, results in a\n", tbeginfit, std::min(ncorr*dt, tendfit));
    }

    tStart = 0;
    if (bPrint)
    {
        printf("COR:%11s%11s%11s%11s%11s%11s%11s\n",
               "Fit from", "Integral", "Tail Value", "Sum (ps)", " a1 (ps)",
               (effnNparams(fitfn) >= 2) ? " a2 ()" : "",
               (effnNparams(fitfn) >= 3) ? " a3 (ps)" : "");
    }

    snew(sig, ncorr);

    if (tbeginfit > 0)
    {
        jmax = 3;
    }
    else
    {
        jmax = 1;
    }
    for (j = 0; ((j < jmax) && (tStart < tendfit) && (tStart < ncorr*dt)); j++)
    {
        /* Estimate the correlation time for better fitting */
        c_start     = -1;
        ct_estimate = 0;
        for (i = 0; (i < ncorr) && (dt*i < tStart || c1[i] > 0); i++)
        {
            if (c_start < 0)
            {
                if (dt*i >= tStart)
                {
                    c_start     = c1[i];
                    ct_estimate = 0.5*c1[i];
                }
            }
            else
            {
                ct_estimate += c1[i];
            }
        }
        if (c_start > 0)
        {
            ct_estimate *= dt/c_start;
        }
        else
        {
            /* The data is strange, so we need to choose somehting */
            ct_estimate = tendfit;
        }
        if (debug)
        {
            fprintf(debug, "tStart %g ct_estimate: %g\n", tStart, ct_estimate);
        }

        if (fitfn == effnEXPEXP)
        {
            fitparm[0] = 0.002*ncorr*dt;
            fitparm[1] = 0.95;
            fitparm[2] = 0.2*ncorr*dt;
        }
        else
        {
            /* Good initial guess, this increases the probability of convergence */
            fitparm[0] = ct_estimate;
            fitparm[1] = 1.0;
            fitparm[2] = 1.0;
        }

        /* Generate more or less appropriate sigma's */
        for (i = 0; i < ncorr; i++)
        {
            sig[i] = sqrt(ct_estimate+dt*i);
        }

        nf_int    = std::min(ncorr, (int)((tStart+1e-4)/dt));
        sum       = print_and_integrate(debug, nf_int, dt, c1, NULL, 1);
        tail_corr = do_lmfit(ncorr, c1, sig, dt, NULL, tStart, tendfit, oenv,
                             bDebugMode(), fitfn, fitparm, 0, NULL);
        sumtot = sum+tail_corr;
        if (fit && ((jmax == 1) || (j == 1)))
        {
            double mfp[3];
            for (i = 0; (i < 3); i++)
            {
                mfp[i] = fitparm[i];
            }
            for (i = 0; (i < ncorr); i++)
            {
                fit[i] = lmcurves[fitfn](i*dt, mfp);
            }
        }
        if (bPrint)
        {
            printf("COR:%11.4e%11.4e%11.4e%11.4e", tStart, sum, tail_corr, sumtot);
            for (i = 0; (i < effnNparams(fitfn)); i++)
            {
                printf(" %11.4e", fitparm[i]);
            }
            printf("\n");
        }
        tStart += tbeginfit;
    }
    sfree(sig);

    return sumtot;
}
