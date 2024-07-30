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
 * Tests for the WholeMoleculeTransform class.
 *
 * \author berk Hess <hess@kth.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "gromacs/mdlib/wholemoleculetransform.h"

#include <cmath>
#include <cstdlib>

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/domdec/ga2la.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

//! The number of atoms for the test molecule
constexpr int c_numAtoms = 4;
//! The number of bonds, we use 2 bonds so we have one unbound atom
constexpr int c_numBonds = 2;

//! Returns a topology with a three atom molecule that is linearly connected
std::unique_ptr<gmx_mtop_t> triAtomMoleculeSystem()
{
    gmx_moltype_t moltype;

    moltype.atoms.nr = c_numAtoms;

    for (int b = 0; b < c_numBonds; b++)
    {
        const int                parameterType = -1; // value is not actually used
        const std::array<int, 2> atomIndices   = { b, b + 1 };
        moltype.ilist[F_CONNBONDS].push_back(parameterType, atomIndices);
    }

    std::unique_ptr<gmx_mtop_t> mtop = std::make_unique<gmx_mtop_t>();

    mtop->moltype.push_back(moltype);

    gmx_molblock_t molblock;
    molblock.type = 0;
    molblock.nmol = 1;
    mtop->molblock.push_back(molblock);
    mtop->natoms = c_numAtoms;

    return mtop;
}


//! Coordinates, broken over PBC
const std::array<RVec, c_numAtoms> c_coords = {
    { { 2.5, 0.5, 1.5 }, { 0.5, 0.5, 1.5 }, { 1.5, 0.5, 1.5 }, { 2.7, 0.5, 1.5 } }
};

//! A box that works with \p c_coords
const matrix box = { { 3, 0, 0 }, { 0, 2, 0 }, { 0, 0, 2 } };

TEST(WholeMoleculeTransform, MakesMoleculesWhole)
{
    const auto mtop = triAtomMoleculeSystem();

    WholeMoleculeTransform wmt(*mtop, PbcType::Xyz, false);

    wmt.updateForAtomPbcJumps(c_coords, box);

    auto wholeCoords = wmt.wholeMoleculeCoordinates(c_coords, box);

    for (int b = 0; b < c_numBonds; b++)
    {
        EXPECT_FLOAT_EQ(std::abs(wholeCoords[b][XX] - wholeCoords[b + 1][XX]), 1);
    }
}

TEST(WholeMoleculeTransform, HandlesReordering)
{
    const auto mtop = triAtomMoleculeSystem();

    WholeMoleculeTransform wmt(*mtop, PbcType::Xyz, true);

    const std::array<int, c_numAtoms> reordering = { 2, 0, 3, 1 };

    gmx_ga2la_t       ga2la(c_numAtoms, c_numAtoms);
    std::vector<RVec> coords;
    for (int a = 0; a < c_numAtoms; a++)
    {
        ga2la.insert(reordering[a], { a, 0 });

        coords.push_back(c_coords[reordering[a]]);
    }

    wmt.updateAtomOrder(reordering, ga2la);

    wmt.updateForAtomPbcJumps(coords, box);

    auto wholeCoords = wmt.wholeMoleculeCoordinates(coords, box);

    for (int b = 0; b < c_numBonds; b++)
    {
        EXPECT_FLOAT_EQ(
                std::abs(wholeCoords[ga2la.find(b)->la][XX] - wholeCoords[ga2la.find(b + 1)->la][XX]), 1);
    }
}

} // namespace
} // namespace test
} // namespace gmx
