/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
/*! \libinternal \file
 * \brief Defines helper functions for mdrun pertaining to free energy calculations.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_mdrunutility
 * \inlibraryapi
 */

#include "gmxpre.h"

#include "freeenergy.h"

#include <memory>
#include <numeric>

#include "gromacs/mdrun/replicaexchange.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"

namespace gmx
{

int computeFepPeriod(const t_inputrec& inputrec, const ReplicaExchangeParameters& replExParams)
{
    if (inputrec.efep == FreeEnergyPerturbationType::No)
    {
        return 0;
    }

    // Set free energy calculation period as the greatest common
    // denominator of nstdhdl, nstcalcenergy, nstexpanded, replica exchange interval,
    // and AWH nstSampleCoord.
    int nstfep = inputrec.fepvals->nstdhdl;
    if (inputrec.nstcalcenergy > 0)
    {
        nstfep = std::gcd(inputrec.nstcalcenergy, nstfep);
    }
    if (inputrec.bExpanded)
    {
        nstfep = std::gcd(inputrec.expandedvals->nstexpanded, nstfep);
    }
    if (replExParams.exchangeInterval > 0)
    {
        nstfep = std::gcd(replExParams.exchangeInterval, nstfep);
    }
    if (inputrec.bDoAwh)
    {
        nstfep = std::gcd(inputrec.awhParams->nstSampleCoord(), nstfep);
    }
    return nstfep;
}

} // namespace gmx
