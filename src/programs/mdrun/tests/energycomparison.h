/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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
using EnergyTermsToCompare = std::unordered_map<std::string, FloatingPointTolerance>;

/*! \internal
 * \brief Function object to compare all energy terms of reference
 * with all matching terms from test within the given tolerances. */
class EnergyComparison
{
public:
    //! Defaults for energy comparisons
    static EnergyTermsToCompare defaultEnergyTermsToCompare();
    //! Constructor
    EnergyComparison(const EnergyTermsToCompare& energyTermsToCompare);
    /*! \brief Return the names of energies that will be compared
     *
     * This function can be used to provide an input for
     * openEnergyFileToReadTerms().
     *
     * \todo This returns a copy of the keys, which is convenient, but
     * inefficient. Alternatively, this could return a view of the keys
     * from a range rather than a container, but there's no implementation
     * of that in C++11 at the moment. */
    std::vector<std::string> getEnergyNames() const;
    /*! \brief Compare \c reference with \c test within \c
     * energyTermsToCompare_
     *
     * Ignore any key found in either \c reference or \c test that is not
     * found in the other. For all keys found in both frames, compare the
     * values with EXPECT_REAL_EQ_TOL and the given tolerance for that
     * key. */
    void operator()(const EnergyFrame& reference, const EnergyFrame& test) const;

    //! Energy terms to match with given tolerances.
    EnergyTermsToCompare energyTermsToCompare_;
};

/*! \brief Check a subset of the energies found in an energy file
 * against reference data.
 *
 * Opens the energy file, loops over all frames, matching the
 * indicated energies against refdata at the given tolerance.
 *
 * \param[in]  energyFilename        The name of an energy file.
 * \param[in]  energyTermsToCompare  Set of energies to match at given tolerances.
 * \param[in]  checker               Root checker for reference data.
 *
 * \todo This is quite similar to the functionality used in PmeTest,
 * and we should consider reducing the duplication.
 */
void checkEnergiesAgainstReferenceData(const std::string&          energyFilename,
                                       const EnergyTermsToCompare& energyTermsToCompare,
                                       TestReferenceChecker*       checker);

} // namespace test
} // namespace gmx

#endif
