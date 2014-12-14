/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "external/lmfit/lmcurve.h"

#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

/*! \brief Number of parameters for each fitting function */
static const int   nfp_ffn[effnNR] = { 0, 1, 2, 3, 2, 5, 7, 9, 4, 3, 6};

const char        *s_ffn[effnNR+2] = {
    NULL, "none", "exp", "aexp", "exp_exp", "vac",
    "exp5", "exp7", "exp9", "erffit", "effnPRES", NULL, NULL
};
/* We don't allow errest as a choice on the command line */

/*! \brief Long description for each fitting function type */
static const char *longs_ffn[effnNR] = {
    "no fit",
    "y = exp(-x/a1)",
    "y = a2 exp(-x/a1)",
    "y = a2 exp(-x/a1) + (1-a2) exp(-x/a3)",
    "y = exp(-v) (cosh(wv) + 1/w sinh(wv)), v = x/(2 a1), w = sqrt(1 - a2)",
    "y = a1 exp(-x/a2) +  a3 exp(-x/a4) + a5",
    "y = a1 exp(-x/a2) +  a3 exp(-x/a4) + a5 exp(-x/a6) + a7",
    "y = a1 exp(-x/a2) +  a3 exp(-x/a4) + a5 exp(-x/a6) + a7 exp(-x/a8) + a9",
    "y = 1/2*(a1+a2) - 1/2*(a1-a2)*erf( (x-a3) /a4^2)",
    "y = a2 *2*a1*((exp(-x/a1)-1)*(a1/x)+1)+(1-a2)*2*a3*((exp(-x/a3)-1)*(a3/x)+1)",
    "y = (1-a1) * exp(-(x/a3)^a4)*cos(x*a2) + a1*exp(-(x/a5)^a6)"
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

    erfarg  = (x-a[2])/(a[3]*a[3]);
    return (a[0]+a[1])/2-(a[0]-a[1])*gmx_erf(erfarg);
}

/*! \brief Compute y = exp(-x/a0) */
static double lmc_exp_one_parm(double x, const double *a)
{
    return exp(-x/a[0]);
}

/*! \brief Compute y = a1 exp(-x/a0) */
static double lmc_exp_two_parm(double x, const double *a)
{
    return a[1]*exp(-x/a[0]);
}

/*! \brief Compute y = a1 exp(-x/a0) + (1-a1) exp(-x/a2) */
static double lmc_exp_3_parm(double x, const double *a)
{
    double e1, e2;

    e1      = exp(-x/a[0]);
    e2      = exp(-x/a[2]);
    return a[1]*e1 + (1-a[1])*e2;
}

/*! \brief Compute y = a0 exp(-x/a1) + a2 exp(-x/a3) + a4 */
static double lmc_exp_5_parm(double x, const double *a)
{
    double e1, e2;

    e1      = exp(-x/a[1]);
    e2      = exp(-x/a[3]);
    return a[0]*e1 + a[2]*e2 + a[4];
}

/*! \brief Compute y = a0 exp(-x/a1) + a2 exp(-x/a3) + a4 exp(-x/a5) + a6 */
static double lmc_exp_7_parm(double x, const double *a)
{
    double e1, e2, e3;

    e1      = exp(-x/a[1]);
    e2      = exp(-x/a[3]);
    e3      = exp(-x/a[5]);
    return a[0]*e1 + a[2]*e2 + a[4]*e3 + a[6];
}

/*! \brief Compute y = a0 exp(-x/a1) + a2 exp(-x/a3) + a4 exp(-x/a5) + a6 exp(-x/a7) + a8 */
static double lmc_exp_9_parm(double x, const double *a)
{
    double e1, e2, e3, e4;

    e1      = exp(-x/a[1]);
    e2      = exp(-x/a[3]);
    e3      = exp(-x/a[5]);
    e4      = exp(-x/a[7]);
    return a[0]*e1 + a[2]*e2 + a[4]*e3 + a[6]*e4 + a[8];
}

/*! \brief Compute y = (1-a1) * exp(-(x/a3)^a4)*cos(x*a2) + a1*exp(-(x/a5)^a6) */
static double effnPRES_fun(double x, const double *a)
{
    double term1, term2, term3;

    if (a[0] != 0)
    {
        term3 = a[0] * exp(-pow((x/a[4]), a[5]));
    }
    else
    {
        term3 = 0;
    }

    term1 = 1-a[0];
    if (term1 != 0)
    {
        term2 = exp(-pow((x/a[2]), a[3])) * cos(x*a[1]);
    }
    else
    {
        term2 = 0;
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
     *    v = x/(2 a1)
     *    w = sqrt(1 - a2)
     *
     *    For tranverse current autocorrelation functions:
     *       a1 = tau
     *       a2 = 4 tau (eta/rho) k^2
     *
     */

    double y, v, det, omega, omega2, em, ec, es;

    v   = x/(2*a[0]);
    det = 1 - a[1];
    em  = exp(-v);
    if (det != 0)
    {
        omega2  = fabs(det);
        omega   = sqrt(omega2);
        if (det > 0)
        {
            ec = em*0.5*(exp(omega*v)+exp(-omega*v));
            es = em*0.5*(exp(omega*v)-exp(-omega*v))/omega;
        }
        else
        {
            ec = em*cos(omega*v);
            es = em*sin(omega*v)/omega;
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
    /*
       e1 = (exp(-x/a1) - 1)
       e2 = (exp(-x/a3) - 1)
       v1=  2*a1 * (e1*a1/x + 1)
       v2 = 2*a3 * (e2*a3/x + 1)
       fun = a2*v1 + (1 - a2) * v2
     */
    double e1, e2, v1, v2;

    if (a[0] != 0)
    {
        e1 = gmx_expm1(-x/a[0]);
    }
    else
    {
        e1 = 0;
    }
    if (a[2] != 0)
    {
        e2 = gmx_expm1(-x/a[2]);
    }
    else
    {
        e2 = 0;
    }

    if (x > 0)
    {
        v1      = 2*a[0]*(e1*a[0]/x + 1);
        v2      = 2*a[2]*(e2*a[2]/x + 1);
        return a[1]*v1 + (1-a[1])*v2;
    }
    else
    {
        return 0;
    }
}

/*! \brief function type for passing to fitting routine */
typedef double (*t_lmcurve)(double x, const double *a);

/*! \brief array of fitting functions corresponding to the pre-defined types */
t_lmcurve lmcurves[effnNR] = {
    lmc_exp_one_parm, lmc_exp_one_parm, lmc_exp_two_parm,
    lmc_exp_3_parm, lmc_vac_2_parm,
    lmc_exp_5_parm,   lmc_exp_7_parm,
    lmc_exp_9_parm, lmc_erffit,  lmc_errest_3_parm, effnPRES_fun
};

double fit_function(int eFitFn, double *parm, double x)
{
    // TODO: check that eFitFn is within bounds
    return lmcurves[eFitFn](x, parm);
}

/*! \brief lmfit_exp supports up to 3 parameter fitting of exponential functions */
static gmx_bool lmfit_exp(int nfit, double x[], double y[],
                          real ftol, int maxiter,
                          double parm[], gmx_bool bVerbose,
                          int eFitFn)
{
    double             chisq, ochisq;
    gmx_bool           bCont;
    int                i, j;
    lm_control_struct *control;
    lm_status_struct  *status;

    if ((eFitFn < 0) || (eFitFn >= effnNR))
    {
        gmx_fatal(FARGS, "fitfn = %d, should be in 0..%d (%s,%d)",
                  effnNR-1, eFitFn, __FILE__, __LINE__);
    }

    snew(control, 1);
    control->ftol       = ftol;
    control->xtol       = ftol;
    control->gtol       = ftol;
    control->epsilon    = 0.1;
    control->stepbound  = 100;
    control->patience   = maxiter;
    control->scale_diag = 1;
    control->msgfile    = stdout;
    control->verbosity  = (bVerbose ? 1 : 0);
    control->n_maxpri   = nfp_ffn[eFitFn];
    control->m_maxpri   = 0;
    snew(status, 1);
    /* Initial params */
    chisq  = 1e12;
    j      = 0;
    if (bVerbose)
    {
        fprintf(stderr, "%4s  %10s  Parameters\n",
                "Step", "chi^2");
    }
    do
    {
        ochisq = chisq;

        lmcurve(nfp_ffn[eFitFn], parm, nfit, x, y,
                lmcurves[eFitFn], control, status);
        chisq = sqr(status->fnorm);
        if (bVerbose)
        {
            int mmm;
            fprintf(stderr, "%4d  %10.5e", j, chisq);
            for (mmm = 0; (mmm < nfp_ffn[eFitFn]); mmm++)
            {
                fprintf(stderr, "  %10.5e", parm[mmm]);
            }
            fprintf(stderr, "\n");
        }
        j++;
        bCont = ((fabs(ochisq - chisq) > fabs(ftol*chisq)) ||
                 ((ochisq == chisq)));
    }
    while (bCont && (j < maxiter));
    if (bVerbose)
    {
        fprintf(stderr, "\n");
    }

    sfree(control);
    sfree(status);

    return TRUE;
}

/*! \brief See description in header file. */
real do_lmfit(int ndata, real c1[], real sig[], real dt, real x0[],
              real begintimefit, real endtimefit, const output_env_t oenv,
              gmx_bool bVerbose, int eFitFn, double fitparms[], int fix)
{
    FILE    *fp;
    char     buf[32];

    int      i, j, nparm, nfitpnts;
    double   integral, ttt;
    double  *parm, *dparm;
    double  *x, *y, *dy;
#ifdef GMX_DOUBLE
    real     ftol    = 1e-16;
#else
    real     ftol    = 1e-8;
#endif
    int      maxiter = 50;

    if (0 != fix)
    {
        fprintf(stderr, "Using fixed parameters in curve fitting is not working anymore\n");
    }
    nparm = nfp_ffn[eFitFn];
    if (debug)
    {
        fprintf(debug, "There are %d points to fit %d vars!\n", ndata, nparm);
        fprintf(debug, "Fit to function %d from %g through %g, dt=%g\n",
                eFitFn, begintimefit, endtimefit, dt);
    }

    snew(x, ndata);
    snew(y, ndata);
    snew(dy, ndata);

    j = 0;
    for (i = 0; i < ndata; i++)
    {
        ttt = x0 ? x0[i] : dt*i;
        if (ttt >= begintimefit && ttt <= endtimefit)
        {
            x[j] = ttt;
            y[j] = c1[i];

            /* mrqmin does not like sig to be zero */
            if (sig[i] < 1.0e-7)
            {
                dy[j] = 1.0e-7;
            }
            else
            {
                dy[j] = sig[i];
            }
            if (debug)
            {
                fprintf(debug, "j= %d, i= %d, x= %g, y= %g, dy= %g\n",
                        j, i, x[j], y[j], dy[j]);
            }
            j++;
        }
    }
    nfitpnts = j;
    integral = 0;
    if (nfitpnts < nparm)
    {
        fprintf(stderr, "Not enough data points for fitting!\n");
    }
    else
    {
        snew(parm, nparm);
        snew(dparm, nparm);
        if (fitparms)
        {
            for (i = 0; (i < nparm); i++)
            {
                parm[i] = fitparms[i];
            }
        }

        if (!lmfit_exp(nfitpnts, x, y, ftol, maxiter, parm, bVerbose, eFitFn))
        {
            fprintf(stderr, "Fit failed!\n");
        }
        else if (nparm <= 3)
        {
            /* Compute the integral from begintimefit */
            if (nparm == 3)
            {
                integral = (parm[0]*myexp(begintimefit, parm[1],  parm[0]) +
                            parm[2]*myexp(begintimefit, 1-parm[1], parm[2]));
            }
            else if (nparm == 2)
            {
                integral = parm[0]*myexp(begintimefit, parm[1],  parm[0]);
            }
            else if (nparm == 1)
            {
                integral = parm[0]*myexp(begintimefit, 1,  parm[0]);
            }
            else
            {
                gmx_fatal(FARGS, "nparm = %d in file %s, line %d",
                          nparm, __FILE__, __LINE__);
            }

            /* Generate THE output */
            if (bVerbose)
            {
                fprintf(stderr, "FIT: # points used in fit is: %d\n", nfitpnts);
                fprintf(stderr, "FIT: %21s%21s%21s\n",
                        "parm0     ", "parm1 (ps)   ", "parm2 (ps)    ");
                fprintf(stderr, "FIT: ------------------------------------------------------------\n");
                fprintf(stderr, "FIT: %8.3g +/- %8.3g%9.4g +/- %8.3g%8.3g +/- %8.3g\n",
                        parm[0], dparm[0], parm[1], dparm[1], parm[2], dparm[2]);
                fprintf(stderr, "FIT: Integral (calc with fitted function) from %g ps to inf. is: %g\n",
                        begintimefit, integral);

                sprintf(buf, "test%d.xvg", nfitpnts);
                fp = xvgropen(buf, "C(t) + Fit to C(t)", "Time (ps)", "C(t)", oenv);
                fprintf(fp, "# parm0 = %g, parm1 = %g, parm2 = %g\n",
                        parm[0], parm[1], parm[2]);
                for (j = 0; (j < nfitpnts); j++)
                {
                    ttt = x0 ? x0[j] : dt*j;
                    fprintf(fp, "%10.5e  %10.5e  %10.5e\n",
                            ttt, c1[j],
                            lmcurves[eFitFn](ttt, parm));
                }
                xvgrclose(fp);
            }
        }
        if (fitparms)
        {
            for (i = 0; (i < nparm); i++)
            {
                fitparms[i] = parm[i];
            }
        }
        sfree(parm);
        sfree(dparm);
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
    nf_int = min(ncorr, (int)(tendfit/dt));
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
        printf("COR:   %s\n", longs_ffn[fitfn]);
        printf("COR: Fit to correlation function from %6.3f ps to %6.3f ps, results in a\n", tbeginfit, min(ncorr*dt, tendfit));
    }

    tStart = 0;
    if (bPrint)
    {
        printf("COR:%11s%11s%11s%11s%11s%11s%11s\n",
               "Fit from", "Integral", "Tail Value", "Sum (ps)", " a1 (ps)",
               (nfp_ffn[fitfn] >= 2) ? " a2 ()" : "",
               (nfp_ffn[fitfn] >= 3) ? " a3 (ps)" : "");
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

        if (fitfn == effnEXP3)
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

        nf_int    = min(ncorr, (int)((tStart+1e-4)/dt));
        sum       = print_and_integrate(debug, nf_int, dt, c1, NULL, 1);
        tail_corr = do_lmfit(ncorr, c1, sig, dt, NULL, tStart, tendfit, oenv,
                             bDebugMode(), fitfn, fitparm, 0);
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
            for (i = 0; (i < nfp_ffn[fitfn]); i++)
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
