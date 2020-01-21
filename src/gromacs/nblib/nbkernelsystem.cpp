/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements nblib kernel system
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "gmxpre.h"

#include <algorithm>

#include "nbkernelsystem.h"

#include "gromacs/math/matrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"

#include "atomtype.h"
#include "simulationstate.h"

namespace nblib
{

NBKernelSystem::NBKernelSystem(SimulationState simState)
{
    const Topology& topology = simState.topology();
    numAtoms                 = topology.numAtoms();

    charges        = topology.getCharges();
    masses         = topology.getMasses();
    excls          = topology.getGmxExclusions();
    atomInfoAllVdw = topology.getAtomInfoAllVdw();

    std::vector<std::tuple<real, real>> nblibNonbonded = topology.getNonbondedParameters();
    nonbondedParameters.resize(numAtoms * numAtoms * 2, 0);
    //! This needs to be filled with all the unique nonbonded params
    //! Todo: Get unique mapping for nonbonded parameters so this can be correctly filled;
    //!       currently this only works this a single nonbonded parameter
    nonbondedParameters[0] = std::get<0>(nblibNonbonded[0]);
    nonbondedParameters[1] = std::get<1>(nblibNonbonded[0]);

    coordinates = simState.coordinates();
    GMX_RELEASE_ASSERT(numAtoms == int(coordinates.size()),
                       "Number of coordinates should match number of atoms in topology");
    velocities = simState.velocities();
    GMX_RELEASE_ASSERT(numAtoms == int(velocities.size()),
                       "Number of velocities should match number of atoms in topology");

    //! Todo: Refactor put_atoms_in_box so that this transformation is not needed
    fillLegacyMatrix(simState.box().matrix(), box);
    put_atoms_in_box(PbcType::Xyz, box, coordinates);

    atomTypes.resize(topology.numAtoms());
    //! This needs to be filled with the atomTypes that correspond to the nonbonded params
    //! Todo: Get unique mapping for atom types so this can be correctly filled;
    //!       currently this only works this a single atom type
    std::fill(std::begin(atomTypes), std::end(atomTypes), 0);
}

} // namespace nblib
