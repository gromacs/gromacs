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
/*! \inpublicapi \file
 * \brief
 * Implements functionality to compute virials from force data
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "virials.h"

#include <algorithm>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/calcvir.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/box.h"
#include "nblib/exception.h"

namespace nblib
{

void computeVirialTensor(gmx::ArrayRef<const Vec3> coordinates,
                         gmx::ArrayRef<const Vec3> forces,
                         gmx::ArrayRef<const Vec3> shiftVectors,
                         gmx::ArrayRef<const Vec3> shiftForces,
                         const Box&                box,
                         gmx::ArrayRef<real>       virialOutput)
{
    if (virialOutput.size() != 9)
    {
        throw InputException("Virial array size incorrect, should be 9");
    }

    // set virial output array to zero
    std::fill(virialOutput.begin(), virialOutput.end(), 0.0);
    // use legacy tensor format
    rvec* virial = reinterpret_cast<rvec*>(virialOutput.data());

    // compute the virials from surrounding boxes
    const rvec* fshift    = as_rvec_array(shiftForces.data());
    const rvec* shift_vec = as_rvec_array(shiftVectors.data());
    calc_vir(numShiftVectors, shift_vec, fshift, virial, false, box.legacyMatrix());

    // calculate partial virial, for local atoms only, based on short range
    // NOTE: GROMACS uses a variable called mdatoms.homenr which is basically number of atoms in current processor
    //       Used numAtoms from the coordinate array in its place
    auto        numAtoms = coordinates.size();
    const rvec* f        = as_rvec_array(forces.data());
    const rvec* x        = as_rvec_array(coordinates.data());
    calc_vir(numAtoms, x, f, virial, false, box.legacyMatrix());
}

} // namespace nblib
