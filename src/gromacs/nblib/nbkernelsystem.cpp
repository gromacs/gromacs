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

#include "nbkernelsystem.h"

#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"

#include "coordinates.h"

namespace nblib
{

NBKernelSystem::NBKernelSystem(const int multiplicationFactor)
{
    nonbondedParameters.resize(numAtoms*numAtoms*2, 0);
    nonbondedParameters[0] = c6Param;
    nonbondedParameters[1] = c12Param;

    generateCoordinates(multiplicationFactor, &coordinates, box);
    put_atoms_in_box(epbcXYZ, box, coordinates);

    int numAtoms = coordinates.size();
    // GMX_RELEASE_ASSERT(numAtoms % numAtomsInMolecule == 0, "Coordinates should match whole molecules");

    // atomTypes.resize(numAtoms);
    // charges.resize(numAtoms);
    // atomInfoAllVdw.resize(numAtoms);
    // snew(excls.index, numAtoms + 1);
    // snew(excls.a, numAtoms*numAtomsInMolecule);
    // excls.index[0] = 0;

    // for (int atom = 0; atom < numAtoms; atom++)
    // {
    //     atomTypes[atom] = atomType;
    //     charges[atom]   = atomCharge;
    //     SET_CGINFO_HAS_VDW(atomInfoAllVdw[atom]);
    //     SET_CGINFO_HAS_Q(atomInfoAllVdw[atom]);

    //     const int firstAtomInMolecule = atom - (atom % numAtomsInMolecule);
    //     for (int atomJ = 0; atomJ < numAtomsInMolecule; atomJ++)
    //     {
    //         excls.a[atom*numAtomsInMolecule + atomJ] = firstAtomInMolecule + atomJ;
    //     }
    //     excls.index[atom + 1] = (atom + 1)*numAtomsInMolecule;
    // }
}

} // namespace nblib
