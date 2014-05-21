/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#include "fit_acf.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/smalloc.h"
#include "expfit.h"

real fit_acf(int ncorr, int fitfn, const output_env_t oenv, gmx_bool bVerbose,
             real tbeginfit, real tendfit, real dt, real c1[], real *fit)
{
    real        fitparm[3];
    real        tStart, tail_corr, sum, sumtot = 0, c_start, ct_estimate, *sig;
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
            for (i = 0; (i < ncorr); i++)
            {
                fit[i] = fit_function(fitfn, fitparm, i*dt);
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
