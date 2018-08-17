/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017,2018, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include <string>

#include "gromacs/utility/stringutil.h"

#include "moduletest.h"

namespace gmx
{
namespace test
{
namespace
{

//! This type holds input integrators.
using CoulombType = const char *;

//! Test fixture parametrized on integrators
class DistributedChargesTest : public gmx::test::MdrunTestFixture,
                               public ::testing::WithParamInterface<CoulombType>
{
};

TEST_P(DistributedChargesTest, Works)
{
    const int         nsteps        = 1;
    const float       timestep      = 0.001;
    auto              coulombtype   = GetParam();
    const std::string coulombtypeName(coulombtype);
    SCOPED_TRACE("Coulombtype " + coulombtypeName);

    int               maxWarningsTolerated = 1;

    const std::string theMdpFile = formatString("integrator              = md\n"
                                                "dt                      = %f\n"
                                                "nsteps                  = %d\n"
                                                "niter                   = 1\n"
                                                "nstcalcenergy           = 1\n"
                                                "cutoff-scheme           = Verlet\n"
                                                "pbc                     = xyz\n"
                                                "coulombtype             = %s\n"
                                                "vdw-type                = Cut-off",
                                                timestep, nsteps, coulombtypeName.c_str());

    runner_.useStringAsMdpFile(theMdpFile);

    const std::string inputFile = "nacl";
    runner_.useTopGroAndNdxFromDatabase(inputFile);
    ::gmx::test::CommandLine Caller;
    Caller.append("grompp");
    Caller.addOption("-maxwarn", maxWarningsTolerated);
    EXPECT_EQ(0, runner_.callGrompp(Caller));

    runner_.edrFileName_ = fileManager_.getTemporaryFilePath(inputFile + ".edr");
    ASSERT_EQ(0, runner_.callMdrun());
}

//! Coulombtypes to test
const CoulombType c_coulombtypesToTest [] = {"Ewald", "PME", "P3M-AD"};

INSTANTIATE_TEST_CASE_P(Checking, DistributedChargesTest, ::testing::ValuesIn(c_coulombtypesToTest));

}  // namespace
}  // namespace test
}  // namespace gmx
