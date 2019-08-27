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
 * \brief Declares function for comparing energy frames.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_ENERGYCOMPARISON_H
#define GMX_PROGRAMS_MDRUN_TESTS_ENERGYCOMPARISON_H

#include <string>
#include <unordered_map>

#include "testutils/testasserts.h"

namespace gmx
{

class EnergyFrame;

namespace test
{

class TestReferenceChecker;

//! Convenience type
using EnergyTolerances = std::unordered_map<std::string, FloatingPointTolerance>;

/*! \brief Compare all fields of reference with all matching fields from test
 *
 * Ignore any key found in either \c reference or \c test that is not
 * found in the other. For all keys found in both frames, compare the
 * values with EXPECT_REAL_EQ_TOL and the given tolerance for that
 * key. */
void compareEnergyFrames(const EnergyFrame      &reference,
                         const EnergyFrame      &test,
                         const EnergyTolerances &tolerances);

/*! \brief Check a subset of the energies found in an energy file
 * against reference data.
 *
 * Opens the energy file, loops over all frames, matching the
 * indicated energies against refdata at the given tolerance.
 *
 * \param[in]  energyFilename   The name of an energy file.
 * \param[in]  energiesToMatch  Set of energies to match at given tolerances.
 * \param[in]  checker          Root checker for reference data.
 *
 * \todo This is quite similar to the functionality used in PmeTest,
 * and we should consider reducing the duplication.
 */
void
checkEnergiesAgainstReferenceData(const std::string      &energyFilename,
                                  const EnergyTolerances &energiesToMatch,
                                  TestReferenceChecker   *checker);

}  // namespace test
}  // namespace gmx

#endif
