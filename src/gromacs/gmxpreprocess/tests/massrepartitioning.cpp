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
#include "gmxpre.h"

#include "gromacs/gmxpreprocess/massrepartitioning.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/fileio/warninp.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

void setupMTop(gmx_mtop_t* mtop, ArrayRef<const real> masses, ArrayRef<const int> bonds)
{
    gmx_moltype_t moltype;

    std::vector<int>& iatoms = moltype.ilist[F_CONNBONDS].iatoms;
    for (Index i = 0; i < bonds.ssize(); i += 2)
    {
        iatoms.push_back(0);
        iatoms.push_back(bonds[i]);
        iatoms.push_back(bonds[i + 1]);
    }

    mtop->moltype.push_back(moltype);

    t_atoms& atoms = mtop->moltype[0].atoms;
    init_atom(&atoms);
    atoms.nr = masses.size();
    snew(atoms.atom, atoms.nr);
    for (Index a = 0; a < ssize(masses); a++)
    {
        atoms.atom[a].m = masses[a];
    }

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;
    mtop->natoms           = atoms.nr;
    mtop->finalize();
}

TEST(MassRepartitioning, ValidCaseWorks)
{
    gmx_mtop_t mtop;

    std::vector<real> masses = { 2, 3, 20 };
    std::vector<int>  bonds  = { 0, 2, 1, 2 };
    setupMTop(&mtop, masses, bonds);

    const real massFactor = 3;

    WarningHandler wi(true, 0);

    // Do the repartitioning
    repartitionAtomMasses(&mtop, false, massFactor, &wi);

    ASSERT_EQ(wi.errorCount() + wi.warningCount(), 0);

    EXPECT_FLOAT_EQ(mtop.moltype[0].atoms.atom[0].m, masses[0] * massFactor);
    EXPECT_FLOAT_EQ(mtop.moltype[0].atoms.atom[1].m, masses[0] * massFactor);
    EXPECT_FLOAT_EQ(mtop.moltype[0].atoms.atom[2].m,
                    masses[2] - (2 * massFactor * masses[0] - masses[0] - masses[1]));
}

TEST(MassRepartitioning, UnboundGivesWarning)
{
    gmx_mtop_t mtop;

    std::vector<real> masses = { 2, 20 };
    std::vector<int>  bonds  = {};
    setupMTop(&mtop, masses, bonds);

    const real massFactor = 3;

    WarningHandler wi(true, 0);

    // Do the repartitioning
    repartitionAtomMasses(&mtop, false, massFactor, &wi);

    EXPECT_EQ(wi.errorCount(), 0);
    EXPECT_EQ(wi.warningCount(), 1);
}

TEST(MassRepartitioning, LightPartnerGivesError)
{
    gmx_mtop_t mtop;

    std::vector<real> masses = { 2, 3, 12 };
    std::vector<int>  bonds  = { 0, 2, 1, 2 };
    setupMTop(&mtop, masses, bonds);

    const real massFactor = 3;

    WarningHandler wi(true, 0);

    // Do the repartitioning
    repartitionAtomMasses(&mtop, false, massFactor, &wi);

    EXPECT_EQ(wi.errorCount(), 1);
    EXPECT_EQ(wi.warningCount(), 0);
}

} // namespace
} // namespace test
} // namespace gmx
