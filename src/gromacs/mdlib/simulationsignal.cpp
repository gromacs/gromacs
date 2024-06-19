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
/*! \internal \file
 *
 * \brief This file defines functions for inter- and intra-simulation
 * signalling by mdrun.
 *
 * This handles details of responding to termination conditions,
 * coordinating checkpoints, and coordinating multi-simulations.
 *
 * \todo Move this to mdrunutility module alongside gathering
 * multi-simulation communication infrastructure there.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdlib
 */

#include "gmxpre.h"

#include "simulationsignal.h"

#include <cstddef>

#include <algorithm>
#include <iterator>

#include "gromacs/gmxlib/network.h"
#include "gromacs/mdrunutility/multisim.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

SimulationSignaller::SimulationSignaller(SimulationSignals*    signals,
                                         const t_commrec*      cr,
                                         const gmx_multisim_t* ms,
                                         bool                  doInterSim,
                                         bool                  doIntraSim) :
    signals_(signals), cr_(cr), ms_(ms), doInterSim_(doInterSim), doIntraSim_(doInterSim || doIntraSim), mpiBuffer_{}
{
}

gmx::ArrayRef<real> SimulationSignaller::getCommunicationBuffer()
{
    if (doIntraSim_)
    {
        std::transform(std::begin(*signals_),
                       std::end(*signals_),
                       std::begin(mpiBuffer_),
                       [](const SimulationSignals::value_type& s) { return s.sig; });

        return mpiBuffer_;
    }
    else
    {
        return {};
    }
}

void SimulationSignaller::signalInterSim()
{
    if (!doInterSim_)
    {
        return;
    }
    // The situations that lead to doInterSim_ == true without a
    // multi-simulation begin active should already have issued an
    // error at mdrun time in release mode, so there's no need for a
    // release-mode assertion.
    GMX_ASSERT(isMultiSim(ms_), "Cannot do inter-simulation signalling without a multi-simulation");
    if (MAIN(cr_))
    {
        // Communicate the signals between the simulations.
        gmx_sum_sim(eglsNR, mpiBuffer_.data(), ms_);
    }
    if (haveDDAtomOrdering(*cr_))
    {
        // Communicate the signals from the main to the others.
        gmx_bcast(eglsNR * sizeof(mpiBuffer_[0]), mpiBuffer_.data(), cr_->mpi_comm_mygroup);
    }
}

void SimulationSignaller::setSignals()
{
    if (!doIntraSim_)
    {
        return;
    }

    SimulationSignals& s = *signals_;
    for (size_t i = 0; i < s.size(); i++)
    {
        if (doInterSim_ || s[i].isLocal)
        {
            // Floating-point reduction of integer values is always
            // exact, so we can use a simple cast.
            signed char gsi = static_cast<signed char>(mpiBuffer_[i]);

            /* Set the communicated signal only when it is non-zero,
             * since a previously set signal might not have been
             * handled immediately. */
            if (gsi != 0)
            {
                s[i].set = gsi;
            }
            // Turn off any local signal now that it has been processed.
            s[i].sig = 0;
        }
    }
}

void SimulationSignaller::finalizeSignals()
{
    signalInterSim();
    setSignals();
}

} // namespace gmx
