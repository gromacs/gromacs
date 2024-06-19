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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "tgroup.h"

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/coupling.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"


real sum_ekin(const t_grpopts* opts, gmx_ekindata_t* ekind, real* dekindlambda, bool bEkinAveVel, bool bScaleEkin)
{
    int           i, j, m, ngtc;
    real          T;
    t_grp_tcstat* tcstat;
    real          nrdf, nd, *ndf;

    ngtc = opts->ngtc;
    ndf  = opts->nrdf;

    T    = 0;
    nrdf = 0;

    clear_mat(ekind->ekin);

    for (i = 0; (i < ngtc); i++)
    {

        nd     = ndf[i];
        tcstat = &ekind->tcstat[i];
        /* Sometimes a group does not have degrees of freedom, e.g.
         * when it consists of shells and virtual sites, then we just
         * set the temperatue to 0 and also neglect the kinetic
         * energy, which should be  zero anyway.
         */

        if (nd > 0)
        {
            if (bEkinAveVel)
            {
                if (!bScaleEkin)
                {
                    /* in this case, kinetic energy is from the current velocities already */
                    msmul(tcstat->ekinf, tcstat->ekinscalef_nhc, tcstat->ekinf);
                }
            }
            else
            {
                /* Calculate the full step Ekin as the average of the half steps */
                for (j = 0; (j < DIM); j++)
                {
                    for (m = 0; (m < DIM); m++)
                    {
                        tcstat->ekinf[j][m] = 0.5
                                              * (tcstat->ekinh[j][m] * tcstat->ekinscaleh_nhc
                                                 + tcstat->ekinh_old[j][m]);
                    }
                }
            }
            m_add(tcstat->ekinf, ekind->ekin, ekind->ekin);

            tcstat->Th = calc_temp(trace(tcstat->ekinh), nd);
            tcstat->T  = calc_temp(trace(tcstat->ekinf), nd);

            /* after the scaling factors have been multiplied in, we can remove them */
            if (bEkinAveVel)
            {
                tcstat->ekinscalef_nhc = 1.0;
            }
            else
            {
                tcstat->ekinscaleh_nhc = 1.0;
            }
        }
        else
        {
            tcstat->T  = 0;
            tcstat->Th = 0;
        }
        T += nd * tcstat->T;
        nrdf += nd;
    }
    if (nrdf > 0)
    {
        T /= nrdf;
    }
    if (dekindlambda)
    {
        if (bEkinAveVel)
        {
            *dekindlambda = ekind->dekindl;
        }
        else
        {
            *dekindlambda = 0.5 * (ekind->dekindl + ekind->dekindl_old);
        }
    }
    return T;
}
