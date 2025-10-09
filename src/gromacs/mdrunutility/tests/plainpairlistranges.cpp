/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Tests PlainPairlistRanges
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "gromacs/mdrunutility/plainpairlistranges.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/topology.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{
namespace
{

std::unique_ptr<gmx_mtop_t> testSystem()
{
    std::unique_ptr<gmx_mtop_t> mtop = std::make_unique<gmx_mtop_t>();

    mtop->moltype.resize(2);
    {
        t_atoms& atoms = mtop->moltype[0].atoms;
        atoms.nr       = 2;
        snew(atoms.atom, atoms.nr);
        atoms.atom[0].m = 14.00670;
        atoms.atom[1].m = 1.00800;
    }
    {
        t_atoms& atoms = mtop->moltype[1].atoms;
        atoms.nr       = 2;
        snew(atoms.atom, atoms.nr);
        atoms.atom[0].m = 15.99940;
        atoms.atom[1].m = 1.00800;
    }
    mtop->molblock.resize(2);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;
    mtop->molblock[1].type = 1;
    mtop->molblock[1].nmol = 5;

    return mtop;
}

TEST(PlainPairlistRanges, RmsdDistance)
{
    const auto mtop = testSystem();
    t_inputrec ir;
    ir.eI                         = IntegrationAlgorithm::MD;
    ir.delta_t                    = 0.002;
    ir.nstlist                    = 25;
    ir.ensembleTemperatureSetting = EnsembleTemperatureSetting::Constant;
    ir.ensembleTemperature        = 298;

    PlainPairlistRanges ppr(*mtop, ir);

    const std::optional<real> rmsdDistance = ppr.rmsdDistance();

    ASSERT_EQ(rmsdDistance.has_value(), true);

    FloatingPointTolerance tolerance(relativeToleranceAsFloatingPoint(1.0, 1e-5));
    EXPECT_REAL_EQ_TOL(rmsdDistance.value(), 0.0776441, tolerance);
}

} // namespace
} // namespace test
} // namespace gmx
