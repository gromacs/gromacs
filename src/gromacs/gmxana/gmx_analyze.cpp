/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <filesystem>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

struct gmx_output_env_t;

/* must correspond to char *avbar_opt[] declared in main() */
enum
{
    avbarSEL,
    avbarNONE,
    avbarSTDDEV,
    avbarERROR,
    avbar90,
    avbarNR
};

static void power_fit(int n, int nset, real** val, real* t)
{
    real *x, *y, quality, a, b, r;
    int   s, i;

    snew(x, n);
    snew(y, n);

    if (t[0] > 0)
    {
        for (i = 0; i < n; i++)
        {
            if (t[0] > 0)
            {
                x[i] = std::log(t[i]);
            }
        }
    }
    else
    {
        fprintf(stdout,
                "First time is not larger than 0, using index number as time for power fit\n");
        for (i = 0; i < n; i++)
        {
            x[i] = std::log1p(static_cast<real>(i));
        }
    }

    for (s = 0; s < nset; s++)
    {
        for (i = 0; i < n && val[s][i] >= 0; i++)
        {
            y[i] = std::log(val[s][i]);
        }
        if (i < n)
        {
            fprintf(stdout, "Will power fit up to point %d, since it is not larger than 0\n", i);
        }
        lsq_y_ax_b(i, x, y, &a, &b, &r, &quality);
        fprintf(stdout, "Power fit set %3d:  error %.3f  a %g  b %g\n", s + 1, quality, a, std::exp(b));
    }

    sfree(y);
    sfree(x);
}

static real cosine_content(int nhp, int n, const real* y)
/* Assumes n equidistant points */
{
    double fac, cosyint, yyint;
    int    i;

    if (n < 2)
    {
        return 0;
    }

    fac = M_PI * nhp / (n - 1);

    cosyint = 0;
    yyint   = 0;
    for (i = 0; i < n; i++)
    {
        cosyint += std::cos(fac * i) * y[i];
        yyint += y[i] * y[i];
    }

    return 2 * cosyint * cosyint / (n * yyint);
}

static void plot_coscont(const char* ccfile, int n, int nset, real** val, const gmx_output_env_t* oenv)
{
    FILE* fp;
    int   s;
    real  cc;

    fp = xvgropen(ccfile, "Cosine content", "set / half periods", "cosine content", oenv);

    for (s = 0; s < nset; s++)
    {
        cc = cosine_content(s + 1, n, val[s]);
        fprintf(fp, " %d %g\n", s + 1, cc);
        fprintf(stdout, "Cosine content of set %d with %.1f periods: %g\n", s + 1, 0.5 * (s + 1), cc);
    }
    fprintf(stdout, "\n");

    xvgrclose(fp);
}

static void regression_analysis(int n, gmx_bool bXYdy, real* x, int nset, real** val)
{
    real S, chi2, a, b, da, db, r = 0;

    if (bXYdy || (nset == 1))
    {
        printf("Fitting data to a function f(x) = ax + b\n");
        printf("Minimizing residual chi2 = Sum_i w_i [f(x_i) - y_i]2\n");
        printf("Error estimates will be given if w_i (sigma) values are given\n");
        printf("(use option -xydy).\n\n");
        if (bXYdy)
        {
            lsq_y_ax_b_error(n, x, val[0], val[1], &a, &b, &da, &db, &r, &S);
        }
        else
        {
            lsq_y_ax_b(n, x, val[0], &a, &b, &r, &S);
        }
        chi2 = gmx::square((n - 2) * S);
        printf("Chi2                    = %g\n", chi2);
        printf("S (Sqrt(Chi2/(n-2))     = %g\n", S);
        printf("Correlation coefficient = %.1f%%\n", 100 * r);
        printf("\n");
        if (bXYdy)
        {
            printf("a    = %g +/- %g\n", a, da);
            printf("b    = %g +/- %g\n", b, db);
        }
        else
        {
            printf("a    = %g\n", a);
            printf("b    = %g\n", b);
        }
    }
    else
    {
        double chi2, *a, **xx, *y;
        int    i, j;

        snew(y, n);
        snew(xx, nset - 1);
        for (j = 0; (j < nset - 1); j++)
        {
            snew(xx[j], n);
        }
        for (i = 0; (i < n); i++)
        {
            y[i] = val[0][i];
            for (j = 1; (j < nset); j++)
            {
                xx[j - 1][i] = val[j][i];
            }
        }
        snew(a, nset - 1);
        chi2 = multi_regression(nullptr, n, y, nset - 1, xx, a);
        printf("Fitting %d data points in %d sets\n", n, nset - 1);
        printf("chi2 = %g\n", chi2);
        printf("A =");
        for (i = 0; (i < nset - 1); i++)
        {
            printf("  %g", a[i]);
            sfree(xx[i]);
        }
        printf("\n");
        sfree(xx);
        sfree(y);
        sfree(a);
    }
}

static void histogram(const char* distfile, real binwidth, int n, int nset, real** val, const gmx_output_env_t* oenv)
{
    FILE*    fp;
    int      i, s;
    double   minval, maxval;
    int      nbin;
    int64_t* histo;

    minval = val[0][0];
    maxval = val[0][0];
    for (s = 0; s < nset; s++)
    {
        for (i = 0; i < n; i++)
        {
            if (val[s][i] < minval)
            {
                minval = val[s][i];
            }
            else if (val[s][i] > maxval)
            {
                maxval = val[s][i];
            }
        }
    }

    minval = binwidth * std::floor(minval / binwidth);
    maxval = binwidth * std::ceil(maxval / binwidth);
    if (minval != 0)
    {
        minval -= binwidth;
    }
    maxval += binwidth;

    nbin = gmx::roundToInt(((maxval - minval) / binwidth) + 1);
    fprintf(stderr, "Making distributions with %d bins\n", nbin);
    snew(histo, nbin);
    fp = xvgropen(distfile, "Distribution", "", "", oenv);
    for (s = 0; s < nset; s++)
    {
        for (i = 0; i < nbin; i++)
        {
            histo[i] = 0;
        }
        for (i = 0; i < n; i++)
        {
            histo[gmx::roundToInt((val[s][i] - minval) / binwidth)]++;
        }
        for (i = 0; i < nbin; i++)
        {
            fprintf(fp, " %g  %g\n", minval + i * binwidth, static_cast<double>(histo[i]) / (n * binwidth));
        }
        if (s < nset - 1)
        {
            fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
        }
    }
    xvgrclose(fp);
}

static int real_comp(const void* a, const void* b)
{
    real dif = *reinterpret_cast<const real*>(a) - *reinterpret_cast<const real*>(b);

    if (dif < 0)
    {
        return -1;
    }
    else if (dif > 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

static void average(const char* avfile, int avbar_opt, int n, int nset, real** val, real* t)
{
    FILE*  fp;
    int    i, s, edge = 0;
    double av, var, err;
    real*  tmp = nullptr;

    fp = gmx_ffopen(avfile, "w");
    if ((avbar_opt == avbarERROR) && (nset == 1))
    {
        avbar_opt = avbarNONE;
    }
    if (avbar_opt != avbarNONE)
    {
        if (avbar_opt == avbar90)
        {
            snew(tmp, nset);
            fprintf(fp, "@TYPE xydydy\n");
            edge = gmx::roundToInt(nset * 0.05);
            fprintf(stdout,
                    "Errorbars: discarding %d points on both sides: %d%%"
                    " interval\n",
                    edge,
                    gmx::roundToInt(100. * (nset - 2 * edge) / nset));
        }
        else
        {
            fprintf(fp, "@TYPE xydy\n");
        }
    }

    for (i = 0; i < n; i++)
    {
        av = 0;
        for (s = 0; s < nset; s++)
        {
            av += val[s][i];
        }
        av /= nset;
        fprintf(fp, " %g %g", t[i], av);
        var = 0;
        if (avbar_opt != avbarNONE)
        {
            if (avbar_opt == avbar90)
            {
                for (s = 0; s < nset; s++)
                {
                    tmp[s] = val[s][i];
                }
                qsort(tmp, nset, sizeof(tmp[0]), real_comp);
                fprintf(fp, " %g %g", tmp[nset - 1 - edge] - av, av - tmp[edge]);
            }
            else
            {
                for (s = 0; s < nset; s++)
                {
                    var += gmx::square(val[s][i] - av);
                }
                if (avbar_opt == avbarSTDDEV)
                {
                    err = std::sqrt(var / nset);
                }
                else
                {
                    err = std::sqrt(var / (nset * (nset - 1)));
                }
                fprintf(fp, " %g", err);
            }
        }
        fprintf(fp, "\n");
    }
    gmx_ffclose(fp);

    if (avbar_opt == avbar90)
    {
        sfree(tmp);
    }
}

/*! \brief Compute final error estimate.
 *
 * See Eqn A17, Hess, JCP 116 (2002) 209-217 for details.
 */
static real optimal_error_estimate(double sigma, const double fitparm[], real tTotal)
{
    // When sigma is zero, the fitparm data can be uninitialized
    if (sigma == 0.0)
    {
        return 0;
    }
    double ss = fitparm[1] * fitparm[0] + (1 - fitparm[1]) * fitparm[2];
    if ((tTotal <= 0) || (ss <= 0))
    {
        fprintf(stderr, "Problem in error estimate: T = %g, ss = %g\n", tTotal, ss);
        return 0;
    }

    return sigma * std::sqrt(2 * ss / tTotal);
}

static void estimate_error(const char*             eefile,
                           int                     nb_min,
                           int                     resol,
                           int                     n,
                           int                     nset,
                           double*                 av,
                           double*                 sig,
                           real**                  val,
                           real                    dt,
                           gmx_bool                bFitAc,
                           gmx_bool                bSingleExpFit,
                           gmx_bool                bAllowNegLTCorr,
                           const gmx_output_env_t* oenv)
{
    FILE*                    fp;
    int                      bs, prev_bs, nbs, nb;
    real                     spacing, nbr;
    int                      s, i, j;
    double                   blav, var;
    std::vector<std::string> leg;
    real *                   tbs, *ybs, rtmp, dens, *fitsig, twooe, tau1_est, tau_sig;
    double                   fitparm[3];
    real                     ee, a, tau1, tau2;

    if (n < 4)
    {
        fprintf(stdout, "The number of points is smaller than 4, can not make an error estimate\n");

        return;
    }

    fp = xvgropen(eefile, "Error estimates", "Block size (time)", "Error estimate", oenv);
    if (output_env_get_print_xvgr_codes(oenv))
    {
        fprintf(fp, "@ subtitle \"using block averaging, total time %g (%d points)\"\n", (n - 1) * dt, n);
    }
    leg.resize(2 * nset);
    xvgrLegend(fp, leg, oenv);

    spacing = std::pow(2.0, 1.0 / resol);
    snew(tbs, n);
    snew(ybs, n);
    snew(fitsig, n);
    for (s = 0; s < nset; s++)
    {
        nbs     = 0;
        prev_bs = 0;
        nbr     = nb_min;
        while (nbr <= n)
        {
            bs = n / static_cast<int>(nbr);
            if (bs != prev_bs)
            {
                nb  = n / bs;
                var = 0;
                for (i = 0; i < nb; i++)
                {
                    blav = 0;
                    for (j = 0; j < bs; j++)
                    {
                        blav += val[s][bs * i + j];
                    }
                    var += gmx::square(av[s] - blav / bs);
                }
                tbs[nbs] = bs * dt;
                if (sig[s] == 0)
                {
                    ybs[nbs] = 0;
                }
                else
                {
                    ybs[nbs] = var / (nb * (nb - 1.0)) * (n * dt) / (sig[s] * sig[s]);
                }
                nbs++;
            }
            nbr *= spacing;
            prev_bs = bs;
        }
        if (sig[s] == 0)
        {
            ee         = 0;
            a          = 1;
            tau1       = 0;
            tau2       = 0;
            fitparm[0] = 0;
            fitparm[1] = 0;
            fitparm[2] = 0;
        }
        else
        {
            for (i = 0; i < nbs / 2; i++)
            {
                rtmp             = tbs[i];
                tbs[i]           = tbs[nbs - 1 - i];
                tbs[nbs - 1 - i] = rtmp;
                rtmp             = ybs[i];
                ybs[i]           = ybs[nbs - 1 - i];
                ybs[nbs - 1 - i] = rtmp;
            }
            /* The initial slope of the normalized ybs^2 is 1.
             * For a single exponential autocorrelation: ybs(tau1) = 2/e tau1
             * From this we take our initial guess for tau1.
             */
            twooe = 2.0 / std::exp(1.0);
            i     = -1;
            do
            {
                i++;
                tau1_est = tbs[i];
            } while (i < nbs - 1 && (ybs[i] > ybs[i + 1] || ybs[i] > twooe * tau1_est));

            if (ybs[0] > ybs[1])
            {
                fprintf(stdout,
                        "Data set %d has strange time correlations:\n"
                        "the std. error using single points is larger than that of blocks of 2 "
                        "points\n"
                        "The error estimate might be inaccurate, check the fit\n",
                        s + 1);
                /* Use the total time as tau for the fitting weights */
                tau_sig = (n - 1) * dt;
            }
            else
            {
                tau_sig = tau1_est;
            }

            if (debug)
            {
                fprintf(debug, "set %d tau1 estimate %f\n", s + 1, tau1_est);
            }

            /* Generate more or less appropriate sigma's,
             * also taking the density of points into account.
             */
            for (i = 0; i < nbs; i++)
            {
                if (i == 0)
                {
                    dens = tbs[1] / tbs[0] - 1;
                }
                else if (i == nbs - 1)
                {
                    dens = tbs[nbs - 1] / tbs[nbs - 2] - 1;
                }
                else
                {
                    dens = 0.5 * (tbs[i + 1] / tbs[i - 1] - 1);
                }
                fitsig[i] = std::sqrt((tau_sig + tbs[i]) / dens);
            }

            if (!bSingleExpFit)
            {
                fitparm[0] = tau1_est;
                fitparm[1] = 0.95;
                /* We set the initial guess for tau2
                 * to halfway between tau1_est and the total time (on log scale).
                 */
                fitparm[2] = std::sqrt(tau1_est * (n - 1) * dt);
                do_lmfit(nbs, ybs, fitsig, 0, tbs, 0, dt * n, oenv, bDebugMode(), effnERREST, fitparm, 0, nullptr);
            }
            if (bSingleExpFit || fitparm[0] < 0 || fitparm[2] < 0 || fitparm[1] < 0
                || (fitparm[1] > 1 && !bAllowNegLTCorr) || fitparm[2] > (n - 1) * dt)
            {
                if (!bSingleExpFit)
                {
                    if (fitparm[2] > (n - 1) * dt)
                    {
                        fprintf(stdout,
                                "Warning: tau2 is longer than the length of the data (%g)\n"
                                "         the statistics might be bad\n",
                                (n - 1) * dt);
                    }
                    else
                    {
                        fprintf(stdout, "a fitted parameter is negative\n");
                    }
                    fprintf(stdout,
                            "invalid fit:  e.e. %g  a %g  tau1 %g  tau2 %g\n",
                            optimal_error_estimate(sig[s], fitparm, n * dt),
                            fitparm[1],
                            fitparm[0],
                            fitparm[2]);
                    /* Do a fit with tau2 fixed at the total time.
                     * One could also choose any other large value for tau2.
                     */
                    fitparm[0] = tau1_est;
                    fitparm[1] = 0.95;
                    fitparm[2] = (n - 1) * dt;
                    fprintf(stdout, "Will fix tau2 at the total time: %g\n", fitparm[2]);
                    do_lmfit(nbs, ybs, fitsig, 0, tbs, 0, dt * n, oenv, bDebugMode(), effnERREST, fitparm, 4, nullptr);
                }
                if (bSingleExpFit || fitparm[0] < 0 || fitparm[1] < 0 || (fitparm[1] > 1 && !bAllowNegLTCorr))
                {
                    if (!bSingleExpFit)
                    {
                        fprintf(stdout, "a fitted parameter is negative\n");
                        fprintf(stdout,
                                "invalid fit:  e.e. %g  a %g  tau1 %g  tau2 %g\n",
                                optimal_error_estimate(sig[s], fitparm, n * dt),
                                fitparm[1],
                                fitparm[0],
                                fitparm[2]);
                    }
                    /* Do a single exponential fit */
                    fprintf(stderr, "Will use a single exponential fit for set %d\n", s + 1);
                    fitparm[0] = tau1_est;
                    fitparm[1] = 1.0;
                    fitparm[2] = 0.0;
                    do_lmfit(nbs, ybs, fitsig, 0, tbs, 0, dt * n, oenv, bDebugMode(), effnERREST, fitparm, 6, nullptr);
                }
            }
            ee   = optimal_error_estimate(sig[s], fitparm, n * dt);
            a    = fitparm[1];
            tau1 = fitparm[0];
            tau2 = fitparm[2];
        }
        fprintf(stdout, "Set %3d:  err.est. %g  a %g  tau1 %g  tau2 %g\n", s + 1, ee, a, tau1, tau2);
        if (output_env_get_xvg_format(oenv) == XvgFormat::Xmgr)
        {
            fprintf(fp, "@ legend string %d \"av %f\"\n", 2 * s, av[s]);
            fprintf(fp,
                    "@ legend string %d \"ee %6g\"\n",
                    2 * s + 1,
                    optimal_error_estimate(sig[s], fitparm, n * dt));
        }
        else if (output_env_get_xvg_format(oenv) == XvgFormat::Xmgrace)
        {
            fprintf(fp, "@ s%d legend \"av %f\"\n", 2 * s, av[s]);
            fprintf(fp, "@ s%d legend \"ee %6g\"\n", 2 * s + 1, optimal_error_estimate(sig[s], fitparm, n * dt));
        }
        for (i = 0; i < nbs; i++)
        {
            fprintf(fp,
                    "%g %g %g\n",
                    tbs[i],
                    sig[s] * std::sqrt(ybs[i] / (n * dt)),
                    sig[s] * std::sqrt(fit_function(effnERREST, fitparm, tbs[i]) / (n * dt)));
        }

        if (bFitAc)
        {
            int    fitlen;
            real * ac, acint;
            double ac_fit[4];

            snew(ac, n);
            for (i = 0; i < n; i++)
            {
                ac[i] = val[s][i] - av[s];
                if (i > 0)
                {
                    fitsig[i] = std::sqrt(static_cast<real>(i));
                }
                else
                {
                    fitsig[i] = 1;
                }
            }
            low_do_autocorr(
                    nullptr, oenv, nullptr, n, 1, -1, &ac, dt, eacNormal, 1, FALSE, TRUE, FALSE, 0, 0, effnNONE);

            fitlen = n / nb_min;

            /* Integrate ACF only up to fitlen/2 to avoid integrating noise */
            acint = 0.5 * ac[0];
            for (i = 1; i <= fitlen / 2; i++)
            {
                acint += ac[i];
            }
            acint *= dt;

            /* Generate more or less appropriate sigma's */
            for (i = 0; i <= fitlen; i++)
            {
                fitsig[i] = std::sqrt(acint + dt * i);
            }

            ac_fit[0] = 0.5 * acint;
            ac_fit[1] = 0.95;
            ac_fit[2] = 10 * acint;
            do_lmfit(n / nb_min, ac, fitsig, dt, nullptr, 0, fitlen * dt, oenv, bDebugMode(), effnEXPEXP, ac_fit, 0, nullptr);
            ac_fit[3] = 1 - ac_fit[1];

            fprintf(stdout,
                    "Set %3d:  ac erest %g  a %g  tau1 %g  tau2 %g\n",
                    s + 1,
                    optimal_error_estimate(sig[s], ac_fit, n * dt),
                    ac_fit[1],
                    ac_fit[0],
                    ac_fit[2]);

            fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            for (i = 0; i < nbs; i++)
            {
                fprintf(fp,
                        "%g %g\n",
                        tbs[i],
                        sig[s] * std::sqrt(fit_function(effnERREST, ac_fit, tbs[i])) / (n * dt));
            }

            sfree(ac);
        }
        if (s < nset - 1)
        {
            fprintf(fp, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
        }
    }
    sfree(fitsig);
    sfree(ybs);
    sfree(tbs);
    xvgrclose(fp);
}

static void luzar_correl(int nn, real* time, int nset, real** val, real temp, gmx_bool bError, real fit_start)
{
    real* kt;
    real  d2;
    int   j;

    please_cite(stdout, "Spoel2006b");

    /* Compute negative derivative k(t) = -dc(t)/dt */
    if (!bError)
    {
        snew(kt, nn);
        compute_derivative(nn, time, val[0], kt);
        for (j = 0; (j < nn); j++)
        {
            kt[j] = -kt[j];
        }
        if (debug)
        {
            d2 = 0;
            for (j = 0; (j < nn); j++)
            {
                d2 += gmx::square(kt[j] - val[3][j]);
            }
            fprintf(debug, "RMS difference in derivatives is %g\n", std::sqrt(d2 / nn));
        }
        analyse_corr(nn, time, val[0], val[2], kt, nullptr, nullptr, nullptr, fit_start, temp);
        sfree(kt);
    }
    else if (nset == 6)
    {
        analyse_corr(nn, time, val[0], val[2], val[4], val[1], val[3], val[5], fit_start, temp);
    }
    else
    {
        printf("Inconsistent input. I need c(t) sigma_c(t) n(t) sigma_n(t) K(t) sigma_K(t)\n");
        printf("Not doing anything. Sorry.\n");
    }
}

static void filter(real flen, int n, int nset, real** val, real dt)
{
    int     f, s, i, j;
    double *filt, sum, vf, fluc, fluctot;

    f = static_cast<int>(flen / (2 * dt));
    snew(filt, f + 1);
    filt[0] = 1;
    sum     = 1;
    for (i = 1; i <= f; i++)
    {
        filt[i] = std::cos(M_PI * dt * i / flen);
        sum += 2 * filt[i];
    }
    for (i = 0; i <= f; i++)
    {
        filt[i] /= sum;
    }
    fprintf(stdout, "Will calculate the fluctuation over %d points\n", n - 2 * f);
    fprintf(stdout, "  using a filter of length %g of %d points\n", flen, 2 * f + 1);
    fluctot = 0;
    for (s = 0; s < nset; s++)
    {
        fluc = 0;
        for (i = f; i < n - f; i++)
        {
            vf = filt[0] * val[s][i];
            for (j = 1; j <= f; j++)
            {
                vf += filt[j] * (val[s][i - f] + val[s][i + f]);
            }
            fluc += gmx::square(val[s][i] - vf);
        }
        fluc /= n - 2 * f;
        fluctot += fluc;
        fprintf(stdout, "Set %3d filtered fluctuation: %12.6e\n", s + 1, std::sqrt(fluc));
    }
    fprintf(stdout, "Overall filtered fluctuation: %12.6e\n", std::sqrt(fluctot / nset));
    fprintf(stdout, "\n");

    sfree(filt);
}

static void do_fit(FILE*                   out,
                   int                     n,
                   gmx_bool                bYdy,
                   int                     ny,
                   real*                   x0,
                   real**                  val,
                   int                     npargs,
                   t_pargs*                ppa,
                   const gmx_output_env_t* oenv,
                   const char*             fn_fitted)
{
    real *  c1 = nullptr, *sig = nullptr;
    double* fitparm;
    real    tendfit, tbeginfit;
    int     i, efitfn, nparm;

    efitfn = get_acffitfn();
    nparm  = effnNparams(efitfn);
    fprintf(out, "Will fit to the following function:\n");
    fprintf(out, "%s\n", effnDescription(efitfn));
    c1 = val[n];
    if (bYdy)
    {
        c1  = val[n];
        sig = val[n + 1];
        fprintf(out, "Using two columns as y and sigma values\n");
    }
    else
    {
        snew(sig, ny);
    }
    if (opt2parg_bSet("-beginfit", npargs, ppa))
    {
        tbeginfit = opt2parg_real("-beginfit", npargs, ppa);
    }
    else
    {
        tbeginfit = x0[0];
    }
    if (opt2parg_bSet("-endfit", npargs, ppa))
    {
        tendfit = opt2parg_real("-endfit", npargs, ppa);
    }
    else
    {
        tendfit = x0[ny - 1];
    }

    snew(fitparm, nparm);
    switch (efitfn)
    {
        case effnEXP1: fitparm[0] = 0.5; break;
        case effnEXP2:
            fitparm[0] = 0.5;
            fitparm[1] = c1[0];
            break;
        case effnEXPEXP:
            fitparm[0] = 1.0;
            fitparm[1] = 0.5 * c1[0];
            fitparm[2] = 10.0;
            break;
        case effnEXP5:
            fitparm[0] = fitparm[2] = 0.5 * c1[0];
            fitparm[1]              = 10;
            fitparm[3]              = 40;
            fitparm[4]              = 0;
            break;
        case effnEXP7:
            fitparm[0] = fitparm[2] = fitparm[4] = 0.33 * c1[0];
            fitparm[1]                           = 1;
            fitparm[3]                           = 10;
            fitparm[5]                           = 100;
            fitparm[6]                           = 0;
            break;
        case effnEXP9:
            fitparm[0] = fitparm[2] = fitparm[4] = fitparm[6] = 0.25 * c1[0];
            fitparm[1]                                        = 0.1;
            fitparm[3]                                        = 1;
            fitparm[5]                                        = 10;
            fitparm[7]                                        = 100;
            fitparm[8]                                        = 0;
            break;
        default:
            fprintf(out, "Warning: don't know how to initialize the parameters\n");
            for (i = 0; (i < nparm); i++)
            {
                fitparm[i] = 1;
            }
    }
    fprintf(out, "Starting parameters:\n");
    for (i = 0; (i < nparm); i++)
    {
        fprintf(out, "a%-2d = %12.5e\n", i + 1, fitparm[i]);
    }
    if (do_lmfit(ny, c1, sig, 0, x0, tbeginfit, tendfit, oenv, bDebugMode(), efitfn, fitparm, 0, fn_fitted) > 0)
    {
        for (i = 0; (i < nparm); i++)
        {
            fprintf(out, "a%-2d = %12.5e\n", i + 1, fitparm[i]);
        }
    }
    else
    {
        fprintf(out, "No solution was found\n");
    }
}

static void print_fitted_function(const char*       fitfile,
                                  const char*       fn_fitted,
                                  gmx_bool          bXYdy,
                                  int               nset,
                                  int               n,
                                  real*             t,
                                  real**            val,
                                  int               npargs,
                                  t_pargs*          ppa,
                                  gmx_output_env_t* oenv)
{
    FILE* out_fit = gmx_ffopen(fitfile, "w");
    if (bXYdy && nset >= 2)
    {
        do_fit(out_fit, 0, TRUE, n, t, val, npargs, ppa, oenv, fn_fitted);
    }
    else
    {
        int s;
        for (s = 0; s < nset; s++)
        {
            std::filesystem::path buf;
            if (fn_fitted)
            {
                buf = gmx::concatenateBeforeExtension(fn_fitted, gmx::formatString("_%d", s));
            }
            do_fit(out_fit, s, FALSE, n, t, val, npargs, ppa, oenv, buf.string().c_str());
        }
    }
    gmx_ffclose(out_fit);
}

int gmx_analyze(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] reads an ASCII file and analyzes data sets.",
        "A line in the input file may start with a time",
        "(see option [TT]-time[tt]) and any number of [IT]y[it]-values may follow.",
        "Multiple sets can also be",
        "read when they are separated by & (option [TT]-n[tt]);",
        "in this case only one [IT]y[it]-value is read from each line.",
        "All lines starting with # and @ are skipped.",
        "All analyses can also be done for the derivative of a set",
        "(option [TT]-d[tt]).[PAR]",

        "All options, except for [TT]-av[tt] and [TT]-power[tt], assume that the",
        "points are equidistant in time.[PAR]",

        "[THISMODULE] always shows the average and standard deviation of each",
        "set, as well as the relative deviation of the third",
        "and fourth cumulant from those of a Gaussian distribution with the same",
        "standard deviation.[PAR]",

        "Option [TT]-ac[tt] produces the autocorrelation function(s).",
        "Be sure that the time interval between data points is",
        "much shorter than the time scale of the autocorrelation.[PAR]",

        "Option [TT]-cc[tt] plots the resemblance of set i with a cosine of",
        "i/2 periods. The formula is::",
        "",
        "  [MATH]2 ([INT][FROM]0[from][TO]T[to][int] y(t) [COS]i [GRK]pi[grk] t[cos] dt)^2 ",
        "  / [INT][FROM]0[from][TO]T[to][int] y^2(t) dt[math]",
        "",
        "This is useful for principal components obtained from covariance",
        "analysis, since the principal components of random diffusion are",
        "pure cosines.[PAR]",

        "Option [TT]-msd[tt] produces the mean square displacement(s).[PAR]",

        "Option [TT]-dist[tt] produces distribution plot(s).[PAR]",

        "Option [TT]-av[tt] produces the average over the sets.",
        "Error bars can be added with the option [TT]-errbar[tt].",
        "The errorbars can represent the standard deviation, the error",
        "(assuming the points are independent) or the interval containing",
        "90% of the points, by discarding 5% of the points at the top and",
        "the bottom.[PAR]",

        "Option [TT]-ee[tt] produces error estimates using block averaging.",
        "A set is divided in a number of blocks and averages are calculated for",
        "each block. The error for the total average is calculated from",
        "the variance between averages of the m blocks B[SUB]i[sub] as follows:",
        "error^2 = [SUM][sum] (B[SUB]i[sub] - [CHEVRON]B[chevron])^2 / (m*(m-1)).",
        "These errors are plotted as a function of the block size.",
        "Also an analytical block average curve is plotted, assuming",
        "that the autocorrelation is a sum of two exponentials.",
        "The analytical curve for the block average is::",
        "",
        "  [MATH]f(t) = [GRK]sigma[grk][TT]*[tt][SQRT]2/T (  [GRK]alpha[grk]   ",
        "  ([GRK]tau[grk][SUB]1[sub] (([EXP]-t/[GRK]tau[grk][SUB]1[sub][exp] - 1) ",
        "  [GRK]tau[grk][SUB]1[sub]/t + 1)) +",
        "  (1-[GRK]alpha[grk]) ([GRK]tau[grk][SUB]2[sub] ",
        "  (([EXP]-t/[GRK]tau[grk][SUB]2[sub][exp] - 1) [GRK]tau[grk][SUB]2[sub]/t + ",
        "  1)))[sqrt][math],",
        "",
        "where T is the total time.",
        "[GRK]alpha[grk], [GRK]tau[grk][SUB]1[sub] and [GRK]tau[grk][SUB]2[sub] are ",
        "obtained by fitting f^2(t) to error^2.",
        "When the actual block average is very close to the analytical curve,",
        "the error is [MATH][GRK]sigma[grk][TT]*[tt][SQRT]2/T (a [GRK]tau[grk][SUB]1[sub] ",
        "+ (1-a) [GRK]tau[grk][SUB]2[sub])[sqrt][math].",
        "The complete derivation is given in",
        "B. Hess, J. Chem. Phys. 116:209-217, 2002.[PAR]",

        "Option [TT]-filter[tt] prints the RMS high-frequency fluctuation",
        "of each set and over all sets with respect to a filtered average.",
        "The filter is proportional to cos([GRK]pi[grk] t/len) where t goes from -len/2",
        "to len/2. len is supplied with the option [TT]-filter[tt].",
        "This filter reduces oscillations with period len/2 and len by a factor",
        "of 0.79 and 0.33 respectively.[PAR]",

        "Option [TT]-g[tt] fits the data to the function given with option",
        "[TT]-fitfn[tt].[PAR]",

        "Option [TT]-power[tt] fits the data to [MATH]b t^a[math], which is accomplished",
        "by fitting to [MATH]a t + b[math] on log-log scale. All points after the first",
        "zero or with a negative value are ignored.[PAR]",

        "Option [TT]-luzar[tt] performs a Luzar & Chandler kinetics analysis",
        "on output from [gmx-hbond]. The input file can be taken directly",
        "from [TT]gmx hbond -ac[tt], and then the same result should be produced.[PAR]",
        "Option [TT]-fitfn[tt] performs curve fitting to a number of different",
        "curves that make sense in the context of molecular dynamics, mainly",
        "exponential curves. More information is in the manual. To check the output",
        "of the fitting procedure the option [TT]-fitted[tt] will print both the",
        "original data and the fitted function to a new data file. The fitting",
        "parameters are stored as comment in the output file."
    };
    static real     tb = -1, te = -1, frac = 0.5, filtlen = 0, binwidth = 0.1, aver_start = 0;
    static gmx_bool bHaveT = TRUE, bDer = FALSE, bSubAv = TRUE, bAverCorr = FALSE, bXYdy = FALSE;
    static gmx_bool bEESEF = FALSE, bEENLC = FALSE, bEeFitAc = FALSE, bPower = FALSE;
    static gmx_bool bIntegrate = FALSE, bRegression = FALSE, bLuzar = FALSE;
    static int      nsets_in = 1, d = 1, nb_min = 4, resol = 10;
    static real     temp = 298.15, fit_start = 1, fit_end = 60;

    /* must correspond to enum avbar* declared at beginning of file */
    static const char* avbar_opt[avbarNR + 1] = { nullptr, "none", "stddev", "error", "90", nullptr };

    t_pargs pa[] = {
        { "-time", FALSE, etBOOL, { &bHaveT }, "Expect a time in the input" },
        { "-b", FALSE, etREAL, { &tb }, "First time to read from set" },
        { "-e", FALSE, etREAL, { &te }, "Last time to read from set" },
        { "-n", FALSE, etINT, { &nsets_in }, "Read this number of sets separated by &" },
        { "-d", FALSE, etBOOL, { &bDer }, "Use the derivative" },
        { "-dp",
          FALSE,
          etINT,
          { &d },
          "HIDDENThe derivative is the difference over this number of points" },
        { "-bw", FALSE, etREAL, { &binwidth }, "Binwidth for the distribution" },
        { "-errbar", FALSE, etENUM, { avbar_opt }, "Error bars for [TT]-av[tt]" },
        { "-integrate",
          FALSE,
          etBOOL,
          { &bIntegrate },
          "Integrate data function(s) numerically using trapezium rule" },
        { "-aver_start", FALSE, etREAL, { &aver_start }, "Start averaging the integral from here" },
        { "-xydy",
          FALSE,
          etBOOL,
          { &bXYdy },
          "Interpret second data set as error in the y values for integrating" },
        { "-regression",
          FALSE,
          etBOOL,
          { &bRegression },
          "Perform a linear regression analysis on the data. If [TT]-xydy[tt] is set a second set "
          "will be interpreted as the error bar in the Y value. Otherwise, if multiple data sets "
          "are present a multilinear regression will be performed yielding the constant A that "
          "minimize [MATH][GRK]chi[grk]^2 = (y - A[SUB]0[sub] x[SUB]0[sub] - A[SUB]1[sub] "
          "x[SUB]1[sub] - ... - A[SUB]N[sub] x[SUB]N[sub])^2[math] where now Y is the first data "
          "set in the input file and x[SUB]i[sub] the others. Do read the information at the "
          "option [TT]-time[tt]." },
        { "-luzar",
          FALSE,
          etBOOL,
          { &bLuzar },
          "Do a Luzar and Chandler analysis on a correlation function and "
          "related as produced by [gmx-hbond]. When in addition the "
          "[TT]-xydy[tt] flag is given the second and fourth column will be "
          "interpreted as errors in c(t) and n(t)." },
        { "-temp",
          FALSE,
          etREAL,
          { &temp },
          "Temperature for the Luzar hydrogen bonding kinetics analysis (K)" },
        { "-fitstart",
          FALSE,
          etREAL,
          { &fit_start },
          "Time (ps) from which to start fitting the correlation functions in order to obtain the "
          "forward and backward rate constants for HB breaking and formation" },
        { "-fitend",
          FALSE,
          etREAL,
          { &fit_end },
          "Time (ps) where to stop fitting the correlation functions in order to obtain the "
          "forward and backward rate constants for HB breaking and formation. Only with "
          "[TT]-gem[tt]" },
        { "-nbmin",
          FALSE,
          etINT,
          { &nb_min },
          "HIDDENMinimum number of blocks for block averaging" },
        { "-resol",
          FALSE,
          etINT,
          { &resol },
          "HIDDENResolution for the block averaging, block size increases with"
          " a factor 2^(1/resol)" },
        { "-eeexpfit",
          FALSE,
          etBOOL,
          { &bEESEF },
          "HIDDENAlways use a single exponential fit for the error estimate" },
        { "-eenlc", FALSE, etBOOL, { &bEENLC }, "HIDDENAllow a negative long-time correlation" },
        { "-eefitac",
          FALSE,
          etBOOL,
          { &bEeFitAc },
          "HIDDENAlso plot analytical block average using a autocorrelation fit" },
        { "-filter",
          FALSE,
          etREAL,
          { &filtlen },
          "Print the high-frequency fluctuation after filtering with a cosine filter of this "
          "length" },
        { "-power", FALSE, etBOOL, { &bPower }, "Fit data to: b t^a" },
        { "-subav", FALSE, etBOOL, { &bSubAv }, "Subtract the average before autocorrelating" },
        { "-oneacf", FALSE, etBOOL, { &bAverCorr }, "Calculate one ACF over all sets" },
    };

    FILE*             out;
    int               n, nlast, s, nset, i, j = 0;
    real **           val, *t, dt, tot, error;
    double *          av, *sig, cum1, cum2, cum3, cum4, db;
    const char *      acfile, *msdfile, *ccfile, *distfile, *avfile, *eefile, *fitfile;
    gmx_output_env_t* oenv;

    t_filenm fnm[] = {
        { efXVG, "-f", "graph", ffREAD },     { efXVG, "-ac", "autocorr", ffOPTWR },
        { efXVG, "-msd", "msd", ffOPTWR },    { efXVG, "-cc", "coscont", ffOPTWR },
        { efXVG, "-dist", "distr", ffOPTWR }, { efXVG, "-av", "average", ffOPTWR },
        { efXVG, "-ee", "errest", ffOPTWR },  { efXVG, "-fitted", "fitted", ffOPTWR },
        { efLOG, "-g", "fitlog", ffOPTWR }
    };
#define NFILE asize(fnm)

    int      npargs;
    t_pargs* ppa;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, npargs, ppa, asize(desc), desc, 0, nullptr, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    acfile   = opt2fn_null("-ac", NFILE, fnm);
    msdfile  = opt2fn_null("-msd", NFILE, fnm);
    ccfile   = opt2fn_null("-cc", NFILE, fnm);
    distfile = opt2fn_null("-dist", NFILE, fnm);
    avfile   = opt2fn_null("-av", NFILE, fnm);
    eefile   = opt2fn_null("-ee", NFILE, fnm);
    if (opt2parg_bSet("-fitfn", npargs, ppa) && acfile == nullptr)
    {
        fitfile = opt2fn("-g", NFILE, fnm);
    }
    else
    {
        fitfile = opt2fn_null("-g", NFILE, fnm);
    }

    val = read_xvg_time(opt2fn("-f", NFILE, fnm),
                        bHaveT,
                        opt2parg_bSet("-b", npargs, ppa),
                        tb,
                        opt2parg_bSet("-e", npargs, ppa),
                        te,
                        nsets_in,
                        &nset,
                        &n,
                        &dt,
                        &t);
    printf("Read %d sets of %d points, dt = %g\n\n", nset, n, dt);

    if (bDer)
    {
        printf("Calculating the derivative as (f[i+%d]-f[i])/(%d*dt)\n\n", d, d);
        n -= d;
        for (s = 0; s < nset; s++)
        {
            for (i = 0; (i < n); i++)
            {
                val[s][i] = (val[s][i + d] - val[s][i]) / (d * dt);
            }
        }
    }

    if (bIntegrate)
    {
        real sum, stddev;

        printf("Calculating the integral using the trapezium rule\n");

        if (bXYdy)
        {
            sum = evaluate_integral(n, t, val[0], val[1], aver_start, &stddev);
            printf("Integral %10.3f +/- %10.5f\n", sum, stddev);
        }
        else
        {
            for (s = 0; s < nset; s++)
            {
                sum = evaluate_integral(n, t, val[s], nullptr, aver_start, &stddev);
                printf("Integral %d  %10.5f  +/- %10.5f\n", s + 1, sum, stddev);
            }
        }
    }

    if (fitfile != nullptr)
    {
        print_fitted_function(
                fitfile, opt2fn_null("-fitted", NFILE, fnm), bXYdy, nset, n, t, val, npargs, ppa, oenv);
    }

    printf("                                      std. dev.    relative deviation of\n");
    printf("                       standard       ---------   cumulants from those of\n");
    printf("set      average       deviation      sqrt(n-1)   a Gaussian distribition\n");
    printf("                                                      cum. 3   cum. 4\n");
    snew(av, nset);
    snew(sig, nset);
    for (s = 0; (s < nset); s++)
    {
        cum1 = 0;
        cum2 = 0;
        cum3 = 0;
        cum4 = 0;
        for (i = 0; (i < n); i++)
        {
            cum1 += val[s][i];
        }
        cum1 /= n;
        for (i = 0; (i < n); i++)
        {
            db = val[s][i] - cum1;
            cum2 += db * db;
            cum3 += db * db * db;
            cum4 += db * db * db * db;
        }
        cum2 /= n;
        cum3 /= n;
        cum4 /= n;
        av[s]  = cum1;
        sig[s] = std::sqrt(cum2);
        if (n > 1)
        {
            error = std::sqrt(cum2 / (n - 1));
        }
        else
        {
            error = 0;
        }
        printf("SS%d  %13.6e   %12.6e   %12.6e      %6.3f   %6.3f\n",
               s + 1,
               av[s],
               sig[s],
               error,
               sig[s] != 0.0 ? cum3 / (sig[s] * sig[s] * sig[s] * std::sqrt(8 / M_PI)) : 0,
               sig[s] != 0.0 ? cum4 / (sig[s] * sig[s] * sig[s] * sig[s] * 3) - 1 : 0);
    }
    printf("\n");

    if (filtlen != 0.0F)
    {
        filter(filtlen, n, nset, val, dt);
    }

    if (msdfile)
    {
        out   = xvgropen(msdfile, "Mean square displacement", "time", "MSD (nm\\S2\\N)", oenv);
        nlast = static_cast<int>(n * frac);
        for (s = 0; s < nset; s++)
        {
            for (j = 0; j <= nlast; j++)
            {
                if (j % 100 == 0)
                {
                    fprintf(stderr, "\r%d", j);
                    fflush(stderr);
                }
                tot = 0;
                for (i = 0; i < n - j; i++)
                {
                    tot += gmx::square(val[s][i] - val[s][i + j]);
                }
                tot /= (n - j);
                fprintf(out, " %g %8g\n", dt * j, tot);
            }
            if (s < nset - 1)
            {
                fprintf(out, "%s\n", output_env_get_print_xvgr_codes(oenv) ? "&" : "");
            }
        }
        xvgrclose(out);
        fprintf(stderr, "\r%d, time=%g\n", j - 1, (j - 1) * dt);
        fflush(stderr);
    }
    if (ccfile)
    {
        plot_coscont(ccfile, n, nset, val, oenv);
    }

    if (distfile)
    {
        histogram(distfile, binwidth, n, nset, val, oenv);
    }
    if (avfile)
    {
        average(avfile, nenum(avbar_opt), n, nset, val, t);
    }
    if (eefile)
    {
        estimate_error(eefile, nb_min, resol, n, nset, av, sig, val, dt, bEeFitAc, bEESEF, bEENLC, oenv);
    }
    if (bPower)
    {
        power_fit(n, nset, val, t);
    }

    if (acfile != nullptr)
    {
        if (bSubAv)
        {
            for (s = 0; s < nset; s++)
            {
                for (i = 0; i < n; i++)
                {
                    val[s][i] -= av[s];
                }
            }
        }
        do_autocorr(acfile, oenv, "Autocorrelation", n, nset, val, dt, eacNormal, bAverCorr);
    }

    if (bRegression)
    {
        regression_analysis(n, bXYdy, t, nset, val);
    }

    if (bLuzar)
    {
        luzar_correl(n, t, nset, val, temp, bXYdy, fit_start);
    }

    view_all(oenv, NFILE, fnm);

    return 0;
}
