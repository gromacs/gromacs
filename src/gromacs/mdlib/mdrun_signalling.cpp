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
/*! \internal \file
 *
 * \brief This file defines functions for inter-rank signalling by mdrun.
 *
 * This handles details of responding to termination conditions,
 * coordinating checkpoints, and coordinating multi-simulations.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "mdrun_signalling.h"

#include <algorithm>

#include "gromacs/legacyheaders/md_support.h"
#include "gromacs/legacyheaders/types/commrec.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

void init_global_signals(struct gmx_signalling_t *gs, const t_commrec *cr,
                         const t_inputrec *ir, int repl_ex_nst)
{
    int i;

    if (MULTISIM(cr))
    {
        gs->nstms = multisim_nstsimsync(cr, ir, repl_ex_nst);
        if (debug)
        {
            fprintf(debug, "Syncing simulations for checkpointing and termination every %d steps\n", gs->nstms);
        }
    }
    else
    {
        gs->nstms = 1;
    }

    for (i = 0; i < eglsNR; i++)
    {
        gs->sig[i] = 0;
        gs->set[i] = 0;
    }
}

gmx::ArrayRef<real>
prepareSignalBuffer(struct gmx_signalling_t *gs)
{
    if (gs)
    {
        gmx::ArrayRef<char> sig(gs->sig);
        gmx::ArrayRef<real> temp(gs->mpiBuffer);

        std::copy(sig.begin(), sig.end(), temp.begin());

        return temp;
    }
    else
    {
        return gmx::EmptyArrayRef();
    }
}

void
handleSignals(struct gmx_signalling_t  *gs,
              const t_commrec          *cr,
              bool                      bInterSimGS)
{
    /* Is the signal in one simulation independent of other simulations? */
    bool bIsSignalLocal[eglsNR] = { false, false, true };

    if (!gs)
    {
        return;
    }

    if (MULTISIM(cr) && bInterSimGS)
    {
        if (MASTER(cr))
        {
            /* Communicate the signals between the simulations */
            gmx_sum_sim(eglsNR, gs->mpiBuffer, cr->ms);
        }
        /* Communicate the signals from the master to the others */
        gmx_bcast(eglsNR*sizeof(gs->mpiBuffer), gs->mpiBuffer, cr);
    }
    for (int i = 0; i < eglsNR; i++)
    {
        if (bInterSimGS || bIsSignalLocal[i])
        {
            /* Set the communicated signal only when it is non-zero,
             * since signals might not be processed at each MD step.
             */
            char gsi = (gs->mpiBuffer[i] >= 0.0 ?
                        static_cast<char>(gs->mpiBuffer[i] + 0.5) :
                        static_cast<char>(gs->mpiBuffer[i] - 0.5));
            if (gsi != 0)
            {
                gs->set[i] = gsi;
            }
            /* Turn off the local signal */
            gs->sig[i] = 0;
        }
    }
}
