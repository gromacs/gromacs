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
 * Implements nblib integrator
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 */
#include "nblib/integrator.h"

#include <cstddef>

#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/particletype.h"
#include "nblib/topology.h"

namespace nblib
{

LeapFrog::LeapFrog(const Topology& topology, const Box& box) : box_(box)
{
    inverseMasses_.resize(topology.numParticles());
    for (int i = 0; i < topology.numParticles(); i++)
    {
        int typeIndex     = topology.getParticleTypeIdOfAllParticles()[i];
        inverseMasses_[i] = 1.0 / topology.getParticleTypes()[typeIndex].mass();
    }
}

LeapFrog::LeapFrog(gmx::ArrayRef<const real> inverseMasses, const Box& box) :
    inverseMasses_(inverseMasses.begin(), inverseMasses.end()), box_(box)
{
}

void LeapFrog::integrate(const real dt, gmx::ArrayRef<Vec3> x, gmx::ArrayRef<Vec3> v, gmx::ArrayRef<const Vec3> f)
{
    for (size_t i = 0; i < x.size(); i++)
    {
        for (int dim = 0; dim < dimSize; dim++)
        {
            v[i][dim] += f[i][dim] * dt * inverseMasses_[i];
            x[i][dim] += v[i][dim] * dt;
        }
    }
}

} // namespace nblib
