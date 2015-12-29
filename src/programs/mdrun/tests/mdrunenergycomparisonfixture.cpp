/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements classes in mdrunenergycomparisonfixture.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "mdrunenergycomparisonfixture.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

#include "energyreader.h"

namespace gmx
{
namespace test
{

MdrunEnergyComparisonFixture::~MdrunEnergyComparisonFixture()
{
}

void MdrunEnergyComparisonFixture::compareFrames(const EnergyFrameInfoPair    &frames,
                                                 std::vector<std::string>      namesOfRequiredFields,
                                                 const FloatingPointTolerance &tolerance)
{
    auto &referenceFrame = frames.first;
    auto &testFrame      = frames.second;

    for (auto const &name : namesOfRequiredFields)
    {
        EXPECT_REAL_EQ_TOL(referenceFrame.getValue(name), testFrame.getValue(name), tolerance)
        << name << " didn't match between reference run " << referenceFrame.getFrameName() << " and test run " << testFrame.getFrameName();
    }
}

namespace
{

//! Database of .mdp strings that supports MdrunEnergyComparisonFixture::prepareSimulation
std::map<const char *, const char *> mdpFileFormatStrings_g
{
    /* Do a highly reproducible simulation on argon
     * - writing frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
     * - with other steps between frame-writing steps
     * - with enough buffer that e.g. a rerun will compute the same potential energy even though it does NS every frame
     * - without constraints
     */
    {
        "argon",
        "rcoulomb = 1.0\n"
        "rvdw = 1.0\n"
        "rlist = -1\n"
        "bd-fric = 1000\n"
        "cutoff-scheme = Verlet\n"
        "verlet-buffer-tolerance = 0.000001\n"
        "nsteps = 20\n"
        "nstcalcenergy = -1\n"
        "nstenergy = 5\n"
        "nstlist = 10\n"
        "nstxout = 5\n"
        "nstvout = 5\n"
        "integrator = %s\n"
        "ld-seed = 234262\n"
        "tcoupl = %s\n"
        "ref-t = 298.0\n"
        "tau-t = 1.0\n"
        "tc-grps = System\n"
        "pcoupl = %s\n"
        "ref-p = 1\n"
        "compressibility = 5e-5\n"
    },
    // As for argon, but this water box tests SETTLE/RATTLE constraints.
    {
        "spc216",
        "rcoulomb = 0.8\n"
        "rvdw = 0.8\n"
        "rlist = -1\n"
        "bd-fric = 1000\n"
        "cutoff-scheme = Verlet\n"
        "verlet-buffer-tolerance = 0.000001\n"
        "nsteps = 20\n"
        "nstcalcenergy = -1\n"
        "nstenergy = 5\n"
        "nstlist = 10\n"
        "nstxout = 5\n"
        "nstvout = 5\n"
        "integrator = %s\n"
        "ld-seed = 234262\n"
        "tcoupl = %s\n"
        "ref-t = 298.0\n"
        "tau-t = 1.0\n"
        "tc-grps = System\n"
        "pcoupl = %s\n"
        "ref-p = 1\n"
        "compressibility = 5e-5\n"
    }
};

}       // namespace

//! Help to prepare simulations (since we can't have a constructor with arguments in a test fixture)
void MdrunEnergyComparisonFixture::prepareSimulation(const char *simulationName,
                                                     const char *integrator,
                                                     const char *tcoupl,
                                                     const char *pcoupl)
{
    runner_.useStringAsMdpFile(formatString(mdpFileFormatStrings_g[simulationName], integrator, tcoupl, pcoupl));
    runner_.useTopGroAndNdxFromDatabase(simulationName);
    EXPECT_EQ(0, runner_.callGrompp());
}

} // namespace test
} // namespace gmx
