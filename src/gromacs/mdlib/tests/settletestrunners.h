/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief Declaration of the SETTLE tests runners.
 *
 * Declares the functions that do the buffer management and apply
 * SETTLE constraints ("test runners").
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H
#define GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H

#include "settletestdata.h"

struct t_pbc;

namespace gmx
{
namespace test
{

/*! \brief Apply SETTLE using CPU version of the algorithm
 *
 * Initializes SETTLE object, applies algorithm, destroys the object. The coordinates, velocities
 * and virial are updated in the testData object.
 *
 * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
 * \param[in]     pbc               Periodic boundary setup.
 * \param[in]     updateVelocities  If the velocities should be updated.
 * \param[in]     calcVirial        If the virial should be computed.
 * \param[in]     testDescription   Brief description that will be printed in case of test failure.
 */
void applySettle(SettleTestData*    testData,
                 t_pbc              pbc,
                 bool               updateVelocities,
                 bool               calcVirial,
                 const std::string& testDescription);

/*! \brief Apply SETTLE using GPU version of the algorithm
 *
 * Initializes SETTLE object, copied data to the GPU, applies algorithm, copies the data back,
 * destroys the object. The coordinates, velocities and virial are updated in the testData object.
 *
 * \param[in,out] testData          An object, containing all the data structures needed by SETTLE.
 * \param[in]     pbc               Periodic boundary setup.
 * \param[in]     updateVelocities  If the velocities should be updated.
 * \param[in]     calcVirial        If the virial should be computed.
 * \param[in]     testDescription   Brief description that will be printed in case of test failure.
 */
void applySettleGpu(SettleTestData*    testData,
                    t_pbc              pbc,
                    bool               updateVelocities,
                    bool               calcVirial,
                    const std::string& testDescription);

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_SETTLETESTRUNNERS_H
