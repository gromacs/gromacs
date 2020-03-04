/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
 * \brief
 * Test for MD with orientation restraints
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class OriresTestFixture : public MdrunTestFixture
{
protected:
    OriresTestFixture();
    ~OriresTestFixture() override;
};

OriresTestFixture::OriresTestFixture() {}

OriresTestFixture::~OriresTestFixture() {}

//! Test fixture for mdrun with orires
typedef gmx::test::OriresTestFixture OriresTest;

/* Check whether the orires function works. */
TEST_F(OriresTest, OriresCanRun)
{
    runner_.useTopGroAndNdxFromDatabase("orires_1lvz");
    const std::string mdpContents = R"(
        dt            = 0.002
        nsteps        = 10
        tcoupl        = Berendsen
        tc-grps       = System
        tau-t         = 0.5
        ref-t         = 300
        constraints   = h-bonds
        cutoff-scheme = Verlet
        orire         = Yes
        orire-fitgrp  = backbone
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine oriresCaller;

    // Do an mdrun with ORIRES enabled
    ASSERT_EQ(0, runner_.callMdrun(oriresCaller));
}

} // namespace test
} // namespace gmx
