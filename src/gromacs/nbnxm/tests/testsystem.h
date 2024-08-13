/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Declares the TestSystem class used for testing NBNxM functionality
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/utility/listoflists.h"

namespace gmx
{
enum class LJCombinationRule : int;

namespace test
{

//! Description of the system used for testing.
struct TestSystem
{
    //! The number of energy groups used by the test system
    static constexpr int sc_numEnergyGroups = 3;

    /*! \brief Constructor
     *
     * Generates test system of a cubic box partially filled with 81 water molecules.
     * It has parts with uncharged molecules, normal SPC/E and part with full LJ.
     *
     * It assigns energy groups in round-robin style based on the largest number of
     * energy groups that might be being tested. This is not general enough to work
     * if we would extend the number of energy-group cases that we test.
     */
    TestSystem(LJCombinationRule ljCombinationRule);

    //! Returns the absolute value of the largest partial charge of the atoms in the system
    static real maxCharge();

    //! Number of different atom types in test system.
    int numAtomTypes;
    //! Storage for parameters for short range interactions.
    std::vector<real> nonbondedParameters;
    //! Storage for atom type parameters.
    std::vector<int> atomTypes;
    //! Storage for atom partial charges.
    std::vector<real> charges;
    //! Atom info
    std::vector<int32_t> atomInfo;
    //! Information about exclusions.
    ListOfLists<int> excls;
    //! Storage for atom positions.
    std::vector<RVec> coordinates;
    //! System simulation box.
    matrix box;
};

} // namespace test

} // namespace gmx
