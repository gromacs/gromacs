/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * \brief
 * Implements ParticleSequencer class
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */

#include "nblib/particlesequencer.h"

#include <cstdio>

#include <algorithm>
#include <stdexcept>

#include "nblib/exception.h"

namespace nblib
{

int ParticleSequencer::operator()(const MoleculeName& moleculeName,
                                  int                 moleculeNr,
                                  const ResidueName&  residueName,
                                  const ParticleName& particleName) const
{
    try
    {
        return data_.at(moleculeName).at(moleculeNr).at(residueName).at(particleName);
    }
    catch (const std::out_of_range& outOfRange)
    {
        if (moleculeName.value() == residueName.value())
        {
            printf("No particle %s in residue %s in molecule %s found\n",
                   particleName.value().c_str(),
                   residueName.value().c_str(),
                   moleculeName.value().c_str());
        }
        else
        {
            printf("No particle %s in molecule %s found\n",
                   particleName.value().c_str(),
                   moleculeName.value().c_str());
        }

        throw InputException(outOfRange.what());
    }
}

void ParticleSequencer::build(const std::vector<std::tuple<Molecule, int>>& moleculesList)
{
    int currentID = 0;
    for (const auto& molNumberTuple : moleculesList)
    {
        const Molecule& molecule = std::get<0>(molNumberTuple);
        const size_t    numMols  = std::get<1>(molNumberTuple);

        auto& moleculeMap = data_[molecule.name()];

        for (size_t i = 0; i < numMols; ++i)
        {
            auto& moleculeNrMap = moleculeMap[i];
            for (int j = 0; j < molecule.numParticlesInMolecule(); ++j)
            {
                moleculeNrMap[molecule.residueName(j)][molecule.particleName(j)] = currentID++;
            }
        }
    }
}

} // namespace nblib
