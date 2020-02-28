/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * Implements nblib integrator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "gmxpre.h"

#include "integrator.h"

#include "gromacs/pbcutil/pbc.h"

namespace nblib
{

LeapFrog::LeapFrog(SimulationState simulationState) : simulationState_(simulationState)
{
    inverseMasses_.resize(simulationState_.topology().numParticles());
    for (int i = 0; i < simulationState_.topology().numParticles(); i++)
    {
        int typeIndex     = simulationState_.topology().getParticleTypeIdOfAllParticles()[i];
        inverseMasses_[i] = 1.0 / simulationState_.topology().getParticleTypes()[typeIndex].mass();
    }
}

void LeapFrog::integrate(const real dt)
{
    std::vector<gmx::RVec>& x = simulationState_.coordinates();
    std::vector<gmx::RVec>& v = simulationState_.velocities();
    std::vector<gmx::RVec>& f = simulationState_.forces();
    for (size_t i = 0; i < x.size(); i++)
    {
        for (int dim = 0; dim < DIM; dim++)
        {
            v[i][dim] += f[i][dim] * dt * inverseMasses_[i];
            x[i][dim] += v[i][dim] * dt;
        }
    }
    matrix box;
    gmx::fillLegacyMatrix(simulationState_.box().matrix(), box);
    put_atoms_in_box(PbcType::Xyz, box, x);
}

} // namespace nblib
