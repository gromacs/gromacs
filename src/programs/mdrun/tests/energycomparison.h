/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief Declares function for comparing energy frames.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#ifndef GMX_PROGRAMS_MDRUN_TESTS_ENERGYCOMPARISON_H
#define GMX_PROGRAMS_MDRUN_TESTS_ENERGYCOMPARISON_H

#include <string>
#include <unordered_map>
#include <vector>

#include "testutils/testasserts.h"

#include "comparison_helpers.h"

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
    EnergyComparison(const EnergyTermsToCompare& energyTermsToCompare, MaxNumFrames maxNumFrames);
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

private:
    //! Energy terms to match with given tolerances.
    EnergyTermsToCompare energyTermsToCompare_;
    //! How many frames should be compared.
    MaxNumFrames maxNumFrames_ = MaxNumFrames::compareAllFrames();
    /*! \brief The number of frames that have been compared until now
     *
     * This field is mutable because the need to update the flag
     * when checking frames is merely an implementation detail,
     * rather than a proper change of internal state triggered
     * by the caller. */
    mutable unsigned int numComparedFrames_ = 0;
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
 * \param[in]  maxNumEnergyFrames    The maximum number of frames to check
 *
 * \todo This is quite similar to the functionality used in PmeTest,
 * and we should consider reducing the duplication.
 */
void checkEnergiesAgainstReferenceData(const std::string&          energyFilename,
                                       const EnergyTermsToCompare& energyTermsToCompare,
                                       TestReferenceChecker*       checker,
                                       MaxNumFrames                maxNumEnergyFrames);

/*!
 * \brief Check a subset of the energies found in an energy file
 * against reference data.
 *
 * Convenience overload using all frames
 *
 * \see checkEnergiesAgainstReferenceData(const std::string&, const EnergyTermsToCompare&, TestReferenceChecker*, MaxNumFrames)
 */
void checkEnergiesAgainstReferenceData(const std::string&          energyFilename,
                                       const EnergyTermsToCompare& energyTermsToCompare,
                                       TestReferenceChecker*       checker);

} // namespace test
} // namespace gmx

#endif
