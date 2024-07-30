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
#include "gmxpre.h"

#include "gromacs/mdlib/updategroupscog.h"

#include <cstdlib>

#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/updategroups.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! Database of 51 water atom input positions (taken from spc216.gro) for use as test inputs.
gmx::RVec positions[] = { { .130, -.041, -.291 },  { .120, -.056, -.192 },  { .044, -.005, -.327 },
                          { -.854, -.406, .477 },  { -.900, -.334, .425 },  { -.858, -.386, .575 },
                          { .351, -.061, .853 },   { .401, -.147, .859 },   { .416, .016, .850 },
                          { -.067, -.796, .873 },  { -.129, -.811, .797 },  { -.119, -.785, .958 },
                          { -.635, -.312, -.356 }, { -.629, -.389, -.292 }, { -.687, -.338, -.436 },
                          { .321, -.919, .242 },   { .403, -.880, .200 },   { .294, -1.001, .193 },
                          { -.404, .735, .728 },   { -.409, .670, .803 },   { -.324, .794, .741 },
                          { .461, -.596, -.135 },  { .411, -.595, -.221 },  { .398, -.614, -.059 },
                          { -.751, -.086, .237 },  { -.811, -.148, .287 },  { -.720, -.130, .152 },
                          { .202, .285, -.364 },   { .122, .345, -.377 },   { .192, .236, -.278 },
                          { -.230, -.485, .081 },  { -.262, -.391, .071 },  { -.306, -.548, .069 },
                          { .464, -.119, .323 },   { .497, -.080, .409 },   { .540, -.126, .258 },
                          { -.462, .107, .426 },   { -.486, .070, .336 },   { -.363, .123, .430 },
                          { .249, -.077, -.621 },  { .306, -.142, -.571 },  { .233, -.110, -.714 },
                          { -.922, -.164, .904 },  { -.842, -.221, .925 },  { -.971, -.204, .827 },
                          { .382, .700, .480 },    { .427, .610, .477 },    { .288, .689, .513 },
                          { .781, .264, -.113 },   { .848, .203, -.070 },   { .708, .283, -.048 } };

TEST(UpdateGroupsCog, ComputesCogs)
{
    const int settleType     = 0;
    const int atomsPerSettle = NRAL(F_SETTLE);
    const int numAtoms       = sizeof(positions) / sizeof(positions[0]);
    const int numMolecules   = gmx::exactDiv(numAtoms, atomsPerSettle);

    // Set up the topology.
    gmx_mtop_t mtop;

    gmx_moltype_t moltype;
    moltype.atoms.nr         = atomsPerSettle;
    std::vector<int>& iatoms = moltype.ilist[F_SETTLE].iatoms;
    iatoms.push_back(settleType);
    iatoms.push_back(0);
    iatoms.push_back(1);
    iatoms.push_back(2);
    mtop.moltype.push_back(moltype);

    mtop.molblock.resize(1);
    mtop.molblock[0].type = 0;
    mtop.molblock[0].nmol = numMolecules;
    mtop.natoms           = numAtoms;
    mtop.finalize();

    // Set up the SETTLE parameters.
    const real dOH = 0.1;
    const real dHH = 0.1633;
    t_iparams  iparams;
    iparams.settle.doh = dOH;
    iparams.settle.dhh = dHH;
    mtop.ffparams.iparams.push_back(iparams);

    // Run the test
    auto result = gmx::makeUpdateGroupingsPerMoleculeType(mtop);
    ASSERT_EQ(std::holds_alternative<std::string>(result), false);
    auto updateGroupingsPerMoleculeType = std::get<std::vector<RangePartitioning>>(result);
    real temperature                    = 300;

    UpdateGroupsCog updateGroupsCog(mtop, updateGroupingsPerMoleculeType, temperature, numAtoms);

    EXPECT_FLOAT_EQ(updateGroupsCog.maxUpdateGroupRadius(), 0.083887339);

    std::vector<int> globalAtomIndices(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        globalAtomIndices[i] = i;
    }

    // Randomize the atom order
    for (int i = 0; i < numAtoms; i++)
    {
        int a1 = std::rand() % numAtoms;
        int a2 = std::rand() % numAtoms;
        if (a1 != a2)
        {
            std::swap(globalAtomIndices[a1], globalAtomIndices[a2]);
            std::swap(positions[a1], positions[a2]);
        }
    }

    updateGroupsCog.addCogs(globalAtomIndices, positions);

    EXPECT_EQ(updateGroupsCog.numCogs(), numMolecules);

    std::vector<gmx::RVec> cogPerAtom(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        cogPerAtom[globalAtomIndices[i]] = updateGroupsCog.cogForAtom(i);
    }
    gmx::test::TestReferenceData    data;
    gmx::test::TestReferenceChecker checker(data.rootChecker());
    checker.checkSequence(cogPerAtom.begin(), cogPerAtom.end(), "cogPerAtom");

    updateGroupsCog.clear();
    EXPECT_EQ(updateGroupsCog.numCogs(), 0);
}

} // namespace
} // namespace test
} // namespace gmx
