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
/*! \internal \file
 * \brief
 * This implements virials tests
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#include "nblib/virials.h"

#include <algorithm>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/arrayref.h"

#include "nblib/box.h"
#include "nblib/nbnxmsetuphelpers.h"

namespace nblib
{
namespace test
{

TEST(VirialsTest, computeVirialTensorWorks)
{
    std::vector<Vec3> coords = { { 0, 1, 2 }, { 2, 3, 4 } };
    std::vector<Vec3> forces = { { 2, 1, 2 }, { 4, 3, 4 } };
    std::vector<Vec3> shiftForces(gmx::c_numShiftVectors, Vec3(0.0, 1.0, 0.0));
    Box               box(1, 2, 3);
    t_forcerec        forcerec;
    updateForcerec(&forcerec, box.legacyMatrix());
    std::vector<Vec3> shiftVectors(gmx::c_numShiftVectors);
    // copy shift vectors from ForceRec
    std::copy(forcerec.shift_vec.begin(), forcerec.shift_vec.end(), shiftVectors.begin());
    std::vector<real> virialTest(9, 0);
    computeVirialTensor(coords, forces, shiftVectors, shiftForces, box, virialTest);
    std::vector<real> virialRef{ -4, -3, -4, -7, -5, -7, -10, -7, -10 };
    EXPECT_EQ(virialRef, virialTest);
}

} // namespace test

} // namespace nblib
