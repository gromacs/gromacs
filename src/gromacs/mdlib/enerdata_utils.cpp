/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "enerdata_utils.h"

#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

gmx_enerdata_t::gmx_enerdata_t(int numEnergyGroups, int numFepLambdas) :
    grpp(numEnergyGroups),
    enerpart_lambda(numFepLambdas == 0 ? 0 : numFepLambdas + 1),
    dhdlLambda(numFepLambdas == 0 ? 0 : numFepLambdas + 1),
    foreign_grpp(numEnergyGroups)
{
}

static real sum_v(int n, gmx::ArrayRef<const real> v)
{
    real t;
    int  i;

    t = 0.0;
    for (i = 0; (i < n); i++)
    {
        t = t + v[i];
    }

    return t;
}

void sum_epot(const gmx_grppairener_t& grpp, real* epot)
{
    int i;

    /* Accumulate energies */
    epot[F_COUL_SR] = sum_v(grpp.nener, grpp.ener[egCOULSR]);
    epot[F_LJ]      = sum_v(grpp.nener, grpp.ener[egLJSR]);
    epot[F_LJ14]    = sum_v(grpp.nener, grpp.ener[egLJ14]);
    epot[F_COUL14]  = sum_v(grpp.nener, grpp.ener[egCOUL14]);

    /* lattice part of LR doesnt belong to any group
     * and has been added earlier
     */
    epot[F_BHAM] = sum_v(grpp.nener, grpp.ener[egBHAMSR]);

    epot[F_EPOT] = 0;
    for (i = 0; (i < F_EPOT); i++)
    {
        if (i != F_DISRESVIOL && i != F_ORIRESDEV)
        {
            epot[F_EPOT] += epot[i];
        }
    }
}

/* Adds computed dV/dlambda contributions to the requested outputs
 *
 * Also adds the dispersion correction dV/dlambda to the VdW term.
 * Note that kinetic energy and constraint contributions are handled in sum_dkdl()
 */
static void sum_dvdl(gmx_enerdata_t* enerd, const t_lambda& fepvals)
{
    // Add dispersion correction to the VdW component
    enerd->dvdl_lin[efptVDW] += enerd->term[F_DVDL_VDW];

    for (size_t i = 0; i < enerd->enerpart_lambda.size(); i++)
    {
        enerd->dhdlLambda[i] += enerd->term[F_DVDL_VDW];
    }

    enerd->term[F_DVDL] = 0.0;
    for (int i = 0; i < efptNR; i++)
    {
        if (fepvals.separate_dvdl[i])
        {
            /* could this be done more readably/compactly? */
            int index;
            switch (i)
            {
                case (efptMASS): index = F_DKDL; break;
                case (efptCOUL): index = F_DVDL_COUL; break;
                case (efptVDW): index = F_DVDL_VDW; break;
                case (efptBONDED): index = F_DVDL_BONDED; break;
                case (efptRESTRAINT): index = F_DVDL_RESTRAINT; break;
                default: index = F_DVDL; break;
            }
            enerd->term[index] = enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvdl-%s[%2d]: %f: non-linear %f + linear %f\n", efpt_names[i], i,
                        enerd->term[index], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
        else
        {
            enerd->term[F_DVDL] += enerd->dvdl_lin[i] + enerd->dvdl_nonlin[i];
            if (debug)
            {
                fprintf(debug, "dvd-%sl[%2d]: %f: non-linear %f + linear %f\n", efpt_names[0], i,
                        enerd->term[F_DVDL], enerd->dvdl_nonlin[i], enerd->dvdl_lin[i]);
            }
        }
    }
}

void accumulatePotentialEnergies(gmx_enerdata_t* enerd, gmx::ArrayRef<const real> lambda, const t_lambda* fepvals)
{
    sum_epot(enerd->grpp, enerd->term);

    if (fepvals)
    {
        /* Note that sum_dvdl() adds the dispersion correction enerd->dvdl_lin[efptVDW],
         * so sum_dvdl() should be called before computing the foreign lambda energy differences.
         */
        sum_dvdl(enerd, *fepvals);

        /* Sum the foreign lambda energy difference contributions.
         * Note that here we only add the potential energy components.
         * The constraint and kinetic energy components are add after integration
         * by sum_dhdl().
         */
        for (int i = 0; i < fepvals->n_lambda; i++)
        {
            /* note we are iterating over fepvals here!
               For the current lam, dlam = 0 automatically,
               so we don't need to add anything to the
               enerd->enerpart_lambda[0] */

            /* we don't need to worry about dvdl_lin contributions to dE at
               current lambda, because the contributions to the current
               lambda are automatically zeroed */

            double& enerpart_lambda = enerd->enerpart_lambda[i + 1];

            for (gmx::index j = 0; j < lambda.ssize(); j++)
            {
                /* Note that this loop is over all dhdl components, not just the separated ones */
                const double dlam = fepvals->all_lambda[j][i] - lambda[j];

                enerpart_lambda += dlam * enerd->dvdl_lin[j];
            }
        }
    }
}

void accumulateKineticLambdaComponents(gmx_enerdata_t*           enerd,
                                       gmx::ArrayRef<const real> lambda,
                                       const t_lambda&           fepvals)
{
    if (fepvals.separate_dvdl[efptBONDED])
    {
        enerd->term[F_DVDL_BONDED] += enerd->term[F_DVDL_CONSTR];
    }
    else
    {
        enerd->term[F_DVDL] += enerd->term[F_DVDL_CONSTR];
    }

    for (int i = 0; i < fepvals.n_lambda; i++)
    {
        /* note we are iterating over fepvals here!
           For the current lam, dlam = 0 automatically,
           so we don't need to add anything to the
           enerd->enerpart_lambda[0] */

        double& enerpart_lambda = enerd->enerpart_lambda[i + 1];

        /* Note that potential energy terms have been added by sum_epot() -> sum_dvdl() */

        /* Constraints can not be evaluated at foreign lambdas, so we add
         * a linear extrapolation. This is an approximation, but usually
         * quite accurate since constraints change little between lambdas.
         */
        const int    lambdaIndex = (fepvals.separate_dvdl[efptBONDED] ? efptBONDED : efptFEP);
        const double dlam        = fepvals.all_lambda[lambdaIndex][i] - lambda[lambdaIndex];
        enerpart_lambda += dlam * enerd->term[F_DVDL_CONSTR];

        if (!fepvals.separate_dvdl[efptMASS])
        {
            const double dlam = fepvals.all_lambda[efptMASS][i] - lambda[efptMASS];
            enerpart_lambda += dlam * enerd->term[F_DKDL];
        }
    }

    /* The constrain contribution is now included in other terms, so clear it */
    enerd->term[F_DVDL_CONSTR] = 0;
}


void reset_foreign_enerdata(gmx_enerdata_t* enerd)
{
    int i, j;

    /* First reset all foreign energy components.  Foreign energies always called on
       neighbor search steps */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->foreign_grpp.ener[i][j] = 0.0;
        }
    }

    /* potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->foreign_term[i] = 0.0;
    }
}

void reset_dvdl_enerdata(gmx_enerdata_t* enerd)
{
    for (int i = 0; i < efptNR; i++)
    {
        enerd->dvdl_lin[i]    = 0.0;
        enerd->dvdl_nonlin[i] = 0.0;
    }
}

void reset_enerdata(gmx_enerdata_t* enerd)
{
    int i, j;

    /* First reset all energy components. */
    for (i = 0; (i < egNR); i++)
    {
        for (j = 0; (j < enerd->grpp.nener); j++)
        {
            enerd->grpp.ener[i][j] = 0.0_real;
        }
    }

    /* Normal potential energy components */
    for (i = 0; (i <= F_EPOT); i++)
    {
        enerd->term[i] = 0.0_real;
    }
    enerd->term[F_PDISPCORR]      = 0.0_real;
    enerd->term[F_DVDL]           = 0.0_real;
    enerd->term[F_DVDL_COUL]      = 0.0_real;
    enerd->term[F_DVDL_VDW]       = 0.0_real;
    enerd->term[F_DVDL_BONDED]    = 0.0_real;
    enerd->term[F_DVDL_RESTRAINT] = 0.0_real;
    enerd->term[F_DKDL]           = 0.0_real;
    std::fill(enerd->enerpart_lambda.begin(), enerd->enerpart_lambda.end(), 0);
    std::fill(enerd->dhdlLambda.begin(), enerd->dhdlLambda.end(), 0);
    /* reset foreign energy data and dvdl - separate functions since they are also called elsewhere */
    reset_foreign_enerdata(enerd);
    reset_dvdl_enerdata(enerd);
}
