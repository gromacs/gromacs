/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
/*! \internal
 * \file
 * \brief
 * Defines a driver routine for lmfit, and a callback for it to use.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "gmx_lmcurve.h"

#include "config.h"

#include <cmath>
#include <cstdio>

#include "gromacs/utility/basedefinitions.h"

#if HAVE_LMFIT
#    include <lmmin.h>
#    include <lmstruct.h>
#endif

#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#if HAVE_LMFIT

typedef struct
{
    const double* t;
    const double* y;
    const double* dy;
    double (*f)(const double t, const double* par);
} lmcurve_data_struct;

//! Callback function used by lmmin
static void lmcurve_evaluate(const double* par, const int m_dat, const void* data, double* fvec, int* info)
{
    const lmcurve_data_struct* D = reinterpret_cast<const lmcurve_data_struct*>(data);
    for (int i = 0; i < m_dat; i++)
    {
        double dy = D->dy[i];
        if (dy == 0)
        {
            dy = 1;
        }
        fvec[i] = (D->y[i] - D->f(D->t[i], par)) / dy;
    }
    *info = 0;
}

//! Calls lmmin with the given data, with callback function \c f.
static void gmx_lmcurve(const int     n_par,
                        double*       par,
                        const int     m_dat,
                        const double* t,
                        const double* y,
                        const double* dy,
                        double (*f)(double t, const double* par),
                        const lm_control_struct* control,
                        lm_status_struct*        status)
{
    lmcurve_data_struct data = { t, y, dy, f };

    lmmin(n_par, par, m_dat, nullptr, &data, lmcurve_evaluate, control, status);
}

#endif

bool lmfit_exp(int          nfit,
               const double x[],
               const double y[],
               const double dy[],
               double       parm[], // NOLINT(readability-non-const-parameter)
               bool         bVerbose,
               int          eFitFn,
               int          nfix)
{
    if ((eFitFn < 0) || (eFitFn >= effnNR))
    {
        fprintf(stderr, "fitfn = %d, should be in the range 0..%d\n", eFitFn, effnNR - 1);
        return false;
    }
#if HAVE_LMFIT
    double            chisq, ochisq;
    gmx_bool          bCont;
    int               j;
    int               maxiter = 100;
    lm_control_struct control;
    lm_status_struct* status;
    int               nparam = effnNparams(eFitFn);
    int               p2;
    gmx_bool          bSkipLast;

    /* Using default control structure for double precision fitting that
     * comes with the lmfit package (i.e. from the include file).
     */
    control           = lm_control_double;
    control.verbosity = (bVerbose ? 1 : 0);
    control.n_maxpri  = 0;
    control.m_maxpri  = 0;

    snew(status, 1);
    /* Initial params */
    chisq = 1e12;
    j     = 0;
    if (bVerbose)
    {
        printf("%4s  %10s  Parameters\n", "Step", "chi^2");
    }
    /* Check whether we have to skip some params */
    if (nfix > 0)
    {
        do
        {
            p2        = 1 << (nparam - 1);
            bSkipLast = ((p2 & nfix) == p2);
            if (bSkipLast)
            {
                nparam--;
                nfix -= p2;
            }
        } while ((nparam > 0) && (bSkipLast));
        if (bVerbose)
        {
            printf("Using %d out of %d parameters\n", nparam, effnNparams(eFitFn));
        }
    }
    do
    {
        ochisq = chisq;
        gmx_lmcurve(nparam, parm, nfit, x, y, dy, lmcurves[eFitFn], &control, status);
        chisq = gmx::square(status->fnorm);
        if (bVerbose)
        {
            printf("status: fnorm = %g, nfev = %d, userbreak = %d\noutcome = %s\n",
                   status->fnorm,
                   status->nfev,
                   status->userbreak,
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
        bCont = (std::fabs(ochisq - chisq) > std::fabs(control.ftol * chisq));
    } while (bCont && (j < maxiter));

    sfree(status);
#else
    gmx_fatal(FARGS,
              "This build of GROMACS was not configured with support "
              "for lmfit, so the requested fitting cannot be performed. "
              "See the install guide for instructions on how to build "
              "GROMACS with lmfit supported.");
    GMX_UNUSED_VALUE(nfit);
    GMX_UNUSED_VALUE(x);
    GMX_UNUSED_VALUE(y);
    GMX_UNUSED_VALUE(dy);
    GMX_UNUSED_VALUE(parm);
    GMX_UNUSED_VALUE(bVerbose);
    GMX_UNUSED_VALUE(eFitFn);
    GMX_UNUSED_VALUE(nfix);
#endif
    return true;
}
