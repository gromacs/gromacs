/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * \brief
 * Implements PlainPairlistRanges
 *
 * \author Berk Hess <hess@kth.se>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */

#include "plainpairlistranges.h"

#include "gromacs/math/units.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

PlainPairlistRanges::PlainPairlistRanges(const gmx_mtop_t& mtop, const t_inputrec& inputrec) :
    mtop_(mtop), inputrec_(inputrec)
{
}

namespace
{

real computeRmsdDistance(const gmx_mtop_t& mtop, const t_inputrec& inputrec)
{
    GMX_RELEASE_ASSERT(EI_DYNAMICS(inputrec.eI), "Need dynamics");
    GMX_RELEASE_ASSERT(inputrec.ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant,
                       "Need a constant ensemble temperature");

    const real lifetime = (inputrec.nstlist - nonbondedMtsFactor(inputrec)) * inputrec.delta_t;

    double sumInvMass = 0;
    int    numAtoms   = 0;
    for (const auto& mb : mtop.molblock)
    {
        const t_atoms& atoms = mtop.moltype[mb.type].atoms;
        for (int i = 0; i < atoms.nr; i++)
        {
            const t_atom& atom = atoms.atom[i];
            const real    minMass =
                    (inputrec.efep == FreeEnergyPerturbationType::No ? atom.m : std::min(atom.m, atom.mB));
            if (minMass > 0)
            {
                sumInvMass += mb.nmol / minMass;
                numAtoms += mb.nmol;
            }
        }
    }

    const real averageInvMass = sumInvMass / std::max(numAtoms, 1);

    return lifetime * std::sqrt(2 * c_boltz * inputrec.ensembleTemperature * averageInvMass);
}

} // namespace


std::optional<real> PlainPairlistRanges::rmsdDistance() const
{
    if (EI_ENERGY_MINIMIZATION(inputrec_.eI))
    {
        return 0;
    }
    else if (EI_DYNAMICS(inputrec_.eI)
             && inputrec_.ensembleTemperatureSetting == EnsembleTemperatureSetting::Constant)
    {
        // Note that we do not precompute this, as most commonly this method is not called
        return computeRmsdDistance(mtop_, inputrec_);
    }

    return std::nullopt;
}

void PlainPairlistRanges::addRange(const real range)
{
    if (range <= 0)
    {
        GMX_THROW(APIError("PlainPairlistRanges.addRange() expects an argument > 0"));
    }

    ranges_.push_back(range);
}

} // namespace gmx
