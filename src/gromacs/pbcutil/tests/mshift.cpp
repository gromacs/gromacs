/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Tests for the graph functionality to make molecules whole that are broken over pbc.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_pbcutil
 */
#include "gmxpre.h"

#include "gromacs/pbcutil/mshift.h"

#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/gmxassert.h"

#include "testutils/testasserts.h"
#include "testutils/testmatchers.h"

namespace gmx
{
namespace test
{
namespace
{

/* TODO: Actually initialize moltype.atoms.atom when this is converted to C++ */

/*! \brief Returns a molectype for generate a graph
 *
 * This moleculetype has the first and last atom not connected to the other atoms
 * so we test optimizations in the shift code that generates the actual graph only
 * for the atom range from the first to the last connected atoms.
 */
gmx_moltype_t moleculeType()
{
    gmx_moltype_t moltype = {};

    moltype.atoms.nr               = 5;
    moltype.ilist[F_CONSTR].iatoms = { 0, 1, 2 };
    moltype.ilist[F_ANGLES].iatoms = { 1, 2, 1, 3 };

    return moltype;
}

//! Box to go with \p coordinates()
constexpr matrix c_box = { { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 } };

/*! \brief Coordinates for \p moleculeType() broken over PBC
 *
 * The middle 3 atoms all need to be shifted with respect to each other
 * to make the molecule whole. The first and last atom are not connected,
 * so their coordinates are irrelevant for this test.
 */
std::vector<RVec> coordinates()
{
    std::vector<RVec> x;

    x.emplace_back(-1, 0, 3);
    x.emplace_back(1.5, 1.5, 4.5);
    x.emplace_back(1.6, 1.5, 1.5);
    x.emplace_back(1.6, -1.3, 1.5);
    x.emplace_back(1, 4, 2);

    return x;
}

/*! \brief Coordinates for \p moleculeType() made whole
 *
 * These coordinates assume the the periodic image for the molecule
 * is chosen the same as the first connected atom.
 */
std::vector<RVec> coordinatesWhole()
{
    std::vector<RVec> x;

    x.emplace_back(-1, 0, 3);
    x.emplace_back(1.5, 1.5, 4.5);
    x.emplace_back(1.6, 1.5, 4.5);
    x.emplace_back(1.6, 1.7, 4.5);
    x.emplace_back(1, 4, 2);

    return x;
}

//! Tests where (un)shifting works to new coordinate buffers
TEST(MShift, shiftsAndUnshifts)
{
    const gmx_moltype_t     molType = moleculeType();
    const std::vector<RVec> x       = coordinates();
    GMX_RELEASE_ASSERT(int(x.size()) == molType.atoms.nr, "coordinates should match moltype");

    t_graph graph = mk_graph_moltype(molType);
    mk_mshift(nullptr, &graph, PbcType::Xyz, c_box, as_rvec_array(x.data()));

    std::vector<RVec> xShifted(molType.atoms.nr);
    shift_x(&graph, c_box, as_rvec_array(x.data()), as_rvec_array(xShifted.data()));
    EXPECT_THAT(coordinatesWhole(), Pointwise(RVecEq(defaultFloatTolerance()), xShifted));

    std::vector<RVec> xUnshifted(molType.atoms.nr);
    unshift_x(&graph, c_box, as_rvec_array(xUnshifted.data()), as_rvec_array(xShifted.data()));
    EXPECT_THAT(x, Pointwise(RVecEq(defaultFloatTolerance()), xUnshifted));
}

//! Tests where (un)shifting works in place
TEST(MShift, shiftsAndUnshiftsSelf)
{
    const gmx_moltype_t molType = moleculeType();
    std::vector<RVec>   x       = coordinates();
    GMX_RELEASE_ASSERT(int(x.size()) == molType.atoms.nr, "coordinates should match moltype");

    t_graph graph = mk_graph_moltype(molType);
    mk_mshift(nullptr, &graph, PbcType::Xyz, c_box, as_rvec_array(x.data()));

    shift_self(&graph, c_box, as_rvec_array(x.data()));
    EXPECT_THAT(coordinatesWhole(), Pointwise(RVecEq(defaultFloatTolerance()), x));

    unshift_self(&graph, c_box, as_rvec_array(x.data()));
    EXPECT_THAT(coordinates(), Pointwise(RVecEq(defaultFloatTolerance()), x));
}

} // namespace
} // namespace test
} // namespace gmx
