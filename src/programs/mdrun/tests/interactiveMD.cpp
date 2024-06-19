/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Tests utilities for interactive molecular dynamics (IMD) setups.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{

class ImdTestFixture : public MdrunTestFixture, public ::testing::WithParamInterface<const char*>
{
protected:
    ImdTestFixture();
    ~ImdTestFixture() override;
};


ImdTestFixture::ImdTestFixture() {}

ImdTestFixture::~ImdTestFixture() {}


//! Test fixture for mdrun with IMD settings
typedef gmx::test::ImdTestFixture ImdTest;

/* This test checks
 * - whether the IMD-group parameter from the .mdp file is understood,
 * - whether mdrun understands the IMD-related command line parameters
     -imdpull, -imdwait, -imdterm,
 * - whether or not GROMACS was compiled with IMD support, that mdrun finishes
     without error when IMD is enabled in the TPR.
 *
 * TODO In future, consider checking that mdrun does not start IMD
 * when it should/can not.
 */
TEST_P(ImdTest, ImdCanRun)
{
    runner_.useTopGroAndNdxFromDatabase("glycine_vacuo");
    const std::string mdpContents = R"(
        dt            = 0.002
        nsteps        = 2
        tcoupl        = v-rescale
        tc-grps       = System
        tau-t         = 0.5
        ref-t         = 300
        cutoff-scheme = Verlet
        IMD-group     = Heavy_Atoms
        integrator    = %s
    )";
    // Interpolate the integrator selection into the .mdp file
    runner_.useStringAsMdpFile(formatString(mdpContents.c_str(), GetParam()));

    EXPECT_EQ(0, runner_.callGrompp());

    ::gmx::test::CommandLine imdCaller;
    imdCaller.addOption("-imdport", 0); // automatically assign a free port
    imdCaller.append("-imdpull");
    imdCaller.append("-noimdwait"); // cannot use -imdwait: then mdrun would not return control ...
    imdCaller.append("-noimdterm");

    // Do an mdrun with IMD enabled
    ASSERT_EQ(0, runner_.callMdrun(imdCaller));
}

// Check a dynamical integrator and an energy minimizer. No need to
// cover the whole space.
INSTANTIATE_TEST_SUITE_P(WithIntegrator, ImdTest, ::testing::Values("md", "steep"));

} // namespace test
} // namespace gmx
